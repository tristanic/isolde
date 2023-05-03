class IsoldeBenchMarker:
    CRYSTAL_BENCHMARKS = {
        'small': ('3io0', '/A:126'),
        'medium': ('6nak', '/A:27'),
        'large': ('8cjh', '/B:153'),
        'huge': ('5zju', '/e:141')
    }
    CRYO_BENCHMARKS = {
        'small': (('7rzq', '24774'), '/C:959'),
        'medium': (('8ehg', '28147'), '/D:170'),
        'large': (('7nhs', '12339'), '/F:533'),
        'huge': (('7oyb', '13112'), '/51:3488'), 
    }

    SIZES = ('small', 'medium', 'large', 'huge')

    ALL_BENCHMARKS = {
        'crystal': CRYSTAL_BENCHMARKS,
        'cryo-EM': CRYO_BENCHMARKS
    }

    def __init__(self, session, max_size='large', max_coord_updates=150, min_coord_updates=10, max_sim_time=300, output_file=None, warning_dialog=True):
        self.session = session
        self.max_coord_updates=max_coord_updates
        self.min_coord_updates=min_coord_updates
        self.max_sim_time=max_sim_time
        from chimerax.core.commands import run
        run(session, 'isolde start')
        from collections import defaultdict
        self.results = defaultdict(dict)
        self.pending_benchmarks = []
        for size in self.SIZES:
            for type_name, benchmarks in self.ALL_BENCHMARKS.items():
                self.pending_benchmarks.append((type_name, size, benchmarks[size]))
            if size==max_size:
                break
        if output_file is not None:
            out = self.output_file = open(output_file, 'wt')
            print(_opengl_info(self.session), file=out)
            print(_system_summary(), file=out)

        else:
            self.output_file = None
        if warning_dialog:
            from .dialog import generic_warning
            generic_warning('ISOLDE will now run a series of representative simulations for benchmarking purposes. To ensure accurate results, please do not interact with ChimeraX until these are complete.')
        self.run_next_benchmark()

    def run_next_benchmark(self, *_):
        from chimerax.core.commands import run
        from chimerax.core.triggerset import DEREGISTER
        if not len(self.pending_benchmarks):
            self.session.logger.info('ISOLDE: All benchmarks completed')
            if self.output_file is not None:
                self.output_file.close()
            return DEREGISTER
        benchmark_type, size, details = self.pending_benchmarks.pop(0)
        self.current_benchmark_type = benchmark_type
        self.current_benchmark_size = size
        if benchmark_type == 'crystal':
            pdb_id, sel_string = details
            self.current_benchmark_name = pdb_id
            sh = run(self.session, f'open {pdb_id} structurefactors true logInfo f')[0]
            xmapset = sh.map_mgr.xmapsets[0]
            static = xmapset.static_xmaps
            if len(static):
                self.session.models.close(static)
            model = sh.structure
        else:
            pdb_id, emdb_id = details[0]
            self.current_benchmark_name = pdb_id
            sel_string = details[1]
            model = run(self.session, f'open {pdb_id} logInfo f')[0]
            v = run(self.session, f'open {emdb_id} from emdb')[0]
            run(self.session, 'ui tool hide "Volume Viewer"')
            run(self.session, f'clipper assoc #{v.id_string} to #{model.id_string}')
        self.session.logger.info(f'Current benchmark model: {pdb_id} ({size} {benchmark_type})')
        self.session.logger.info('='*25)
        run(self.session, f'addh #{model.id_string} metaldist 1')
        run(self.session, f'isolde select #{model.id_string}')
        br = self._current_runner = BenchMarkRunner(self.session, model, benchmark_type, sel_string, self.max_coord_updates, self.min_coord_updates, self.max_sim_time)
        br.triggers.add_handler('finished', self.benchmark_finished_cb)
        return DEREGISTER



    def benchmark_finished_cb(self, trigger_name, result):
        self.results[self.current_benchmark_type][self.current_benchmark_size] = (self.current_benchmark_name, result)
        if self.output_file is not None:
            out = self.output_file
            print(f'PDB ID:\t{self.current_benchmark_name}', file=out)
            print('='*20, file=out)
            for sel_string, rdata in result.items():
                print(f'Selection string:\t{sel_string}', file=out)
                for val_label, val in result[sel_string].items():
                    print(f'{val_label}:\t{val}', file=out)
                print('-'*10,file=out)
        from chimerax.core.commands import run
        run(self.session, 'close')
        self.session.triggers.add_handler('new frame', self.run_next_benchmark)
    
    


class BenchMarkRunner:
    def __init__(self, session, model, benchmark_type, local_sel_string, max_coord_updates=150, min_coord_updates=10, max_time=300):
        self.session = session
        self.model = model
        self.benchmark_type = benchmark_type
        self.sel_strings = [f'#{model.id_string}', f'#{model.id_string}'+local_sel_string]
        self.max_coord_updates = max_coord_updates
        self.min_coord_updates = min_coord_updates
        self.update_count = 0
        self.max_time = max_time
        from chimerax.core.triggerset import TriggerSet
        t = self.triggers = TriggerSet()
        t.add_trigger('finished')
        from collections import defaultdict
        self.results = defaultdict(dict)
        self.start_benchmark()

    
    def start_benchmark(self):
        sel_string = self.sel_strings.pop(0)
        self.session.logger.info(f'Simulation selection string: {sel_string}')
        self.session.logger.info('-'*25)
        self.current_results = self.results[sel_string]
        from time import time
        self.start_time = self.last_time = time()
        from chimerax.core.commands import run
        self.minimizing=True
        self.first_start=True
        self.update_count = 0
        self.coord_update_times = []
        run(self.session, 'show ~HC')
        run(self.session, f'view {sel_string}')
        run(self.session, f'isolde sim start {sel_string}; sel clear')
        sh = self.session.isolde.sim_handler
        sh.triggers.add_handler('coord update', self._coord_update_cb)
        sh.triggers.add_handler('sim terminated', self._sim_end_cb)
        sim_atom_count = len(self.session.isolde.sim_manager.sim_construct.all_atoms)
        self.current_results['Simulated atom count'] = sim_atom_count
        self.session.logger.info(f'Simulated atom count: {sim_atom_count}')
        platform = sh._context.getPlatform().getName()
        self.current_results['Platform'] = platform
        self.session.logger.info(f'Simulation platform: {platform}')
        if self.benchmark_type=='crystal':
            self.map_tracker = MapUpdateTracker(self.model.parent.map_mgr.xmapsets[0])
        self.graphics_tracker = GraphicsPerformanceTracker(self.session)
    
    def _sim_end_cb(self, *_):
        import numpy
        cut = numpy.array(self.coord_update_times)
        self.current_results['Time per coord update (mean)'] = cut.mean()
        self.current_results['Time per coord update (std)'] = cut.std()
        self.session.logger.info(f'Time per coord update (mean): {cut.mean():.4f} seconds')
        self.session.logger.info(f'Time per coord update (std): {cut.std():.4f} seconds')
        if self.benchmark_type == 'crystal':
            map_times = numpy.array(self.map_tracker.times)
            self.current_results['Time per x-ray map recalculation (mean)'] = map_times.mean()
            self.current_results['Time per x-ray map recalculation (std)'] = map_times.std()
            self.session.logger.info(f'Time per x-ray map recalculation (mean): {map_times.mean():.4f} seconds')
            self.session.logger.info(f'Time per x-ray map recalculation (std): {map_times.std():.4f} seconds')
            self.map_tracker.handler.remove()
        frame_times = numpy.array(self.graphics_tracker.times)
        self.current_results['Time per graphics update (mean)'] = frame_times.mean()
        self.current_results['Time per graphics update (std)'] = frame_times.std()
        self.current_results['Time per graphics update (slowest)'] = frame_times.max()
        self.session.logger.info(f'Time per graphics update (mean): {frame_times.mean()*1000:.1f} ms')
        self.session.logger.info(f'Time per graphics update (std): {frame_times.std()*1000:.1f} ms')
        self.session.logger.info(f'Time per graphics update (slowest): {frame_times.max()*1000:.1f} ms')
        self.session.logger.info('-'*25)
        self.graphics_tracker.handler.remove()
        if not len(self.sel_strings):
            self.triggers.activate_trigger('finished', dict(self.results))
        else:
            self.start_benchmark()

    
    def _coord_update_cb(self, *_):
        self.update_count += 1
        sh = self.session.isolde.sim_handler
        from time import time
        current_time = time()
        if self.first_start:
            self.current_results['Time to first coord update'] = current_time - self.start_time
            self.session.logger.info(f'Time to first coord update: {current_time-self.start_time:.4f} seconds')
            self.first_start = False
            self.last_time = current_time
        elif sh.minimize: 
            pass
        else:
            if self.minimizing:
                self.current_results['Minimization time'] = current_time - self.last_time
                self.session.logger.info(f'Minimization time: {current_time - self.last_time:.4f} seconds')
                self.minimizing = False
            else:
                self.coord_update_times.append(current_time - self.last_time)
            self.last_time = current_time
        if (current_time - self.start_time > self.max_time and self.update_count >= self.min_coord_updates) or self.update_count > self.max_coord_updates:
            from chimerax.core.commands import run
            run(self.session, 'isolde sim stop discard start')

class MapUpdateTracker:
    def __init__(self, xmapset):
        self.last_time = None
        self.times = []
        self.handler = xmapset.triggers.add_handler('maps recalculated', self._recalc_cb)
    
    def _recalc_cb(self, *_):
        from time import time
        current_time = time()
        if self.last_time is not None:
            self.times.append(current_time-self.last_time)
        self.last_time = current_time

class GraphicsPerformanceTracker:
    def __init__(self, session, discard_first=3):
        self._startup_counter = 0
        self._discard_first_n = discard_first
        self.times = []
        self.handler = session.triggers.add_handler('new frame', self._new_frame_cb)
        self._last_time = None
    
    def _new_frame_cb(self, *_):
        self._startup_counter += 1
        if self._startup_counter < self._discard_first_n:
            return
        from time import time
        if self._last_time is None:
            self._last_time = time()
        else:
            current_time = time()
            self.times.append(current_time-self._last_time)
            self._last_time = current_time

def _system_summary():
    import sys
    from chimerax.bug_reporter.bug_reporter_gui import _darwin_info, _win32_info, _linux_info
    if sys.platform == 'win32':
        return _win32_info()
    if sys.platform == 'linux':
        return _linux_info()
    if sys.platform == 'darwin':
        return _darwin_info()

def _opengl_info(session):
    r = session.main_view.render
    try:
        r.make_current()
        lines = ['OpenGL version: ' + r.opengl_version(),
                    'OpenGL renderer: ' + r.opengl_renderer(),
                    'OpenGL vendor: ' + r.opengl_vendor()]
        r.done_current()
    except Exception:
        lines = ['OpenGL version: unknown',
                    'Could not make opengl context current']
    return '\n'.join(lines)


def run_benchmarks(session, max_size='large', output_file=None, warning_dialog=True, max_coord_updates=150, min_coord_updates=10, max_sim_time=300):
    if output_file is None:
        output_file = 'isolde_benchmark.log'
    def _run(*_):
        IsoldeBenchMarker(session, max_size, max_coord_updates, min_coord_updates, max_sim_time, output_file, warning_dialog)
        from chimerax.core.triggerset import DEREGISTER
        return DEREGISTER
    if not hasattr(session, 'isolde'):
        from .cmd import isolde_start
        isolde_start(session, show_splash=False)
        session.isolde.triggers.add_handler('gui started', _run)
    else:
        _run()


def register_isolde_benchmark(logger):
    from chimerax.core.commands import (
        register, CmdDesc,
        EnumOf, PositiveIntArg, PositiveFloatArg, FileNameArg, BoolArg
    )
    desc = CmdDesc(
        synopsis=('Run a series of representative crystallographic and cryo-EM models in ISOLDE to benchmark system performance'),
        keyword=[
            ('max_size', EnumOf(IsoldeBenchMarker.SIZES)),
            ('output_file', FileNameArg),
            ('warning_dialog', BoolArg),
            ('max_coord_updates', PositiveIntArg),
            ('min_coord_updates', PositiveIntArg),
            ('max_sim_time', PositiveFloatArg)
        ]
    )
    register('isolde benchmark', desc, run_benchmarks, logger=logger)




    



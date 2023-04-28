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

    def __init__(self, session, max_size='large'):
        self.session = session
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
        self.run_next_benchmark()

    def run_next_benchmark(self):
        from chimerax.core.commands import run
        run(self.session, 'close')
        if not len(self.pending_benchmarks):
            self.session.logger.info('ISOLDE: All benchmarks completed')
            return
        benchmark_type, size, details = self.pending_benchmarks.pop(0)
        self.current_benchmark_type = benchmark_type
        self.current_benchmark_size = size
        if benchmark_type == 'crystal':
            pdb_id, sel_string = details
            self.current_benchmark_name = pdb_id
            sh = run(self.session, f'open {pdb_id} structurefactors true')[0]
            xmapset = sh.map_mgr.xmapsets[0]
            static = xmapset.static_xmaps
            if len(static):
                self.session.models.close(static)
            model = sh.structure
        else:
            pdb_id, emdb_id = details[0]
            self.current_benchmark_name = pdb_id
            sel_string = details[1]
            model = run(self.session, f'open {pdb_id}')[0]
            v = run(self.session, f'open {emdb_id} from emdb')[0]
            run(self.session, f'clipper assoc #{v.id_string} to #{model.id_string}')
        run(self.session, f'addh #{model.id_string} metaldist 1')
        run(self.session, f'isolde select #{model.id_string}')
        br = self._current_runner = BenchMarkRunner(self.session, model, benchmark_type, sel_string)
        br.triggers.add_handler('finished', self.benchmark_finished_cb)



    def benchmark_finished_cb(self, trigger_name, result):
        self.results[self.current_benchmark_type][self.current_benchmark_size] = (self.current_benchmark_name, result)
        self.run_next_benchmark()

    


class BenchMarkRunner:
    def __init__(self, session, model, benchmark_type, local_sel_string, max_coord_updates=100, max_time=30):
        self.session = session
        self.model = model
        self.benchmark_type = benchmark_type
        self.sel_strings = [f'#{model.id_string}', f'#{model.id_string}'+local_sel_string]
        self.max_coord_updates = max_coord_updates
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
        self.current_results = self.results[sel_string]
        from time import time
        self.start_time = self.last_time = time()
        from chimerax.core.commands import run
        self.minimizing=True
        self.first_start=True
        self.update_count = 0
        self.coord_update_times = []
        run(self.session, 'show ~HC')
        run(self.session, f'view {sel_string}&protein')
        run(self.session, f'isolde sim start {sel_string}; sel clear')
        sh = self.session.isolde.sim_handler
        sh.triggers.add_handler('coord update', self._coord_update_cb)
        sh.triggers.add_handler('sim terminated', self._sim_end_cb)
        if self.benchmark_type=='crystal':
            self.map_tracker = MapUpdateTracker(self.model.parent.map_mgr.xmapsets[0])
    
    def _sim_end_cb(self, *_):
        import numpy
        cut = numpy.array(self.coord_update_times)
        self.current_results['Time per coord update (mean)'] = cut.mean()
        self.current_results['Time per coord update (std)'] = cut.std()
        if self.benchmark_type == 'crystal':
            map_times = numpy.array(self.map_tracker.times)
            self.current_results['Time per x-ray map recalculation (mean)'] = map_times.mean()
            self.current_results['Time per x-ray map recalculation (std)'] = map_times.std()
            self.map_tracker.handler.remove()
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
            self.first_start = False
            self.last_time = current_time
        elif sh.minimize: 
            pass
        else:
            if self.minimizing:
                self.current_results['Minimization time'] = current_time - self.last_time
                self.minimizing = False
            else:
                self.coord_update_times.append(current_time - self.last_time)
            self.last_time = current_time
        if current_time - self.start_time > self.max_time or self.update_count > self.max_coord_updates:
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


        




    



# @Author: Tristan Croll <tic20>
# @Date:   01-Aug-2020
# @Email:  tic20@cam.ac.uk
# @Last modified by:   tic20
# @Last modified time: 03-Aug-2020
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll

global _menu_prepared
_menu_prepared = False

def prepare_isolde_menu(session):
    global _menu_prepared
    if _menu_prepared:
        return
    try:
        main_win = session.ui.main_window

        # from chimerax.core.commands import runscript

        import glob, os
        import importlib.util

        base_dir = os.path.dirname(os.path.abspath(__file__))
        entries = glob.glob(os.path.join(base_dir, '**/*.py'), recursive=True)
        for entry in entries:
            e_path = os.path.split(os.path.relpath(entry, start=base_dir))
            file = e_path[-1]
            if file in ('__init__.py', 'menu.py'):
                continue
            full_path = os.path.abspath(entry)
            try:
                spec = importlib.util.spec_from_file_location(os.path.splitext(file)[0], full_path)
                foo = importlib.util.module_from_spec(spec)

                spec.loader.exec_module(foo)
                tooltip = getattr(foo, 'tooltip', '')
                menu_labels = os.path.normpath(e_path[0]).split(os.sep)
                formatted_labels = ['ISOLDE']
                for label in menu_labels:
                    if len(label):
                        formatted_labels.append(format_label(label))
                menu_entry = format_label(file)
                # def f(session=session):
                #     print('Running script')
                #     foo.run_script(session)
                main_win.add_menu_entry(formatted_labels, menu_entry, item_runner_factory(session, foo), tool_tip=tooltip)
            except Exception as e:
                warn_str = ('Initialisation of ISOLDE menu item "{}" failed: {}').format(
                    '/'.format(menu_labels)+'/'+menu_entry, str(e)
                )
        _menu_prepared=True
    except Exception as e:
        session.logger.warning('Initialisation of ISOLDE menu failed: {}'.format(str(e)))

def item_runner_factory(session, module):
    def f(session=session, module=module):
        module.run_script(session)
    return f

def format_label(label):
    import os
    label = os.path.splitext(label)[0]
    return ' '.join(label.title().split('_'))

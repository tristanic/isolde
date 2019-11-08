_available_tutorials = {
    'Intro to cryo-EM model building': ['intro','cryo_intro','cryo_intro.html'],
    'Intro to crystallographic model building': ['intro', 'crystal_intro', 'crystal_intro.html'],
    'Flexible fitting with adaptive distance restraints': ['bulk_fitting', 'bulk_fitting.html'],
}

def populate_tutorial_combo_box(combo_box):
    import os
    combo_box.blockSignals(True)
    combo_box.clear()
    combo_box.addItem('Tutorials', None)
    base_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
        'docs', 'user', 'tutorials')
    for key, tut_path in _available_tutorials.items():
        combo_box.addItem(key, os.path.join(base_path, *tut_path))
    combo_box.setCurrentIndex(0)
    combo_box.blockSignals(False)

def show_tutorial(session, tut_path):
    from chimerax.help_viewer import show_url
    import pathlib
    show_url(session, pathlib.Path(tut_path).as_uri())

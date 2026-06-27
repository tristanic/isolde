# @Author: Tristan Croll <tic20>
# @Date:   11-Jun-2019
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 11-Jun-2019
# @License: Free for non-commercial use (see license.pdf)
# @Copyright: 2016-2019 Tristan Croll



'''
Dialog boxes for use by the ISOLDE gui
'''

# Reference to the ISOLDE splash screen while it is on display. The splash is
# created with Qt.WindowStaysOnTopHint, so any modal dialog raised while it is
# up would otherwise appear *behind* it (and the render loop that fades it out
# is blocked by the dialog's nested event loop). We therefore dismiss the splash
# before showing any blocking dialog. tool.py registers/clears it here.
_active_splash = None

def register_splash(splash):
    '''
    Record (or clear, with ``None``) the ISOLDE splash screen currently on
    display, so that any blocking dialog can dismiss it first.
    '''
    global _active_splash
    _active_splash = splash

def _dismiss_splash():
    global _active_splash
    splash = _active_splash
    _active_splash = None
    if splash is not None:
        try:
            splash.close()
        except RuntimeError:
            # Underlying Qt object has already been deleted; nothing to do.
            pass

def generic_warning(message):
    _dismiss_splash()
    from Qt.QtWidgets import QMessageBox
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText(message)
    msg.setStandardButtons(QMessageBox.Ok)
    msg.exec()

def choice_warning(message, allow_dont_ask_again=False, yesno=False):
    '''
    Pop up a warning dialog box with the given message, and return True
    if the user wants to go ahead.
    '''
    _dismiss_splash()
    from Qt.QtWidgets import QMessageBox
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msg.setText(message)
    if yesno:
        msg.setStandardButtons(QMessageBox.Yes|QMessageBox.No)
    else:
        msg.setStandardButtons(QMessageBox.Ok|QMessageBox.Cancel)
    if allow_dont_ask_again:
        dag_button = msg.addButton("OK (don't ask again)", QMessageBox.ButtonRole.YesRole)
    reply = msg.exec()
    if allow_dont_ask_again:
        if msg.clickedButton()==dag_button:
            return True, True
        elif yesno:
            if reply==QMessageBox.Yes:
                return True, False
            return False, False
        else:
            if reply==QMessageBox.Ok:
                return True, False
            return False, False
    if yesno:
        if reply == QMessageBox.Yes:
            return True
        return False
    if reply == QMessageBox.Ok:
        return True
    return False

def failed_template_warning(residue):
    '''
    Warning dialog handling the case where a template is not recognised by
    OpenMM when attempting to start a simulation.
    '''
    _dismiss_splash()
    from Qt.QtWidgets import QMessageBox, QPushButton
    msg = QMessageBox()
    msg.setIcon(QMessageBox.Warning)
    msgtext = 'Residue {} {} of chain {} (shown) does not match any template'\
        + ' in the molecular dynamics database. It may be missing atoms (have'\
        + ' you added hydrogens?) or be an unusual residue that has not been'\
        + ' parameterised. Choose what you wish to do with it from the options'\
        + ' below.'
    msg.setText(msgtext.format(residue.name, residue.number, residue.chain_id))

    addh = QPushButton('Add hydrogens and retry')
    msg.addButton(addh, QMessageBox.AcceptRole)
    exclude = QPushButton('Exclude residue from simulations and retry')
    msg.addButton(exclude, QMessageBox.RejectRole)
    abort = QPushButton('Abort')
    msg.addButton(abort, QMessageBox.NoRole)
    msg.exec()
    btn = msg.clickedButton()
    # print("Button: {}".format(btn))
    if btn == addh:
        return "addh"
    if btn == exclude:
        return "exclude"
    return "abort"

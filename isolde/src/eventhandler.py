# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


# vim: set expandtab shiftwidth=4 softtabstop=4:

# Generic wrapper for handling of events. It would also make sense to put
# custom mouse modes here
class EventHandler():
    def __init__(self, session):
        self.session = session
        self.registered_handlers = {}
    # Types of event available in ChimeraX that are useful for ISOLDE
    event_keys = [
        'selection changed', #User has added or removed to an atom selection
        'new frame',
        'shape changed',
        'add models',
        'remove models',
        'frame drawn',
        'graphics update'
        ]

    

    def add_event_handler(self, name, event_type, callback):
        if event_type not in self.event_keys:
            raise Exception('Tried to add handler for unrecognised event')
        if name in self.registered_handlers:
            raise Exception('Handler must have a unique name')
        
        handler = self.session.triggers.add_handler(event_type, callback)
        self.registered_handlers[name] = handler
        return handler
        
        
    def remove_event_handler(self,name):
        if name not in self.registered_handlers:
            raise Exception('Tried to remove unrecognised handler')
        handler = self.registered_handlers[name]
        self.session.triggers.remove_handler(handler)
        self.registered_handlers.pop(name,'Handler not found')
    
    def list_event_handlers(self):
        return self.registered_handlers.keys()
    
    
        

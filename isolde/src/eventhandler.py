# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


# vim: set expandtab shiftwidth=4 softtabstop=4:

# Generic wrapper for handling of events. It would also make sense to put
# custom mouse modes here
class EventHandler():
    def __init__(self, owner):
        self.owner = owner
        self.registered_handlers = {}

    

    def add_event_handler(self, name, event_type, callback):
        
        handler = self.owner.triggers.add_handler(event_type, callback)
        self.registered_handlers[name] = handler
        return handler
        
        
    def remove_event_handler(self,name):
        if name not in self.registered_handlers:
            raise Exception('Tried to remove unrecognised handler')
        handler = self.registered_handlers[name]
        self.owner.triggers.remove_handler(handler)
        self.registered_handlers.pop(name)
    
    def remove_all_handlers(self):
        names = list(self.list_event_handlers())
        for name in names:
            handler = self.registered_handlers[name]
            self.owner.triggers.remove_handler(handler)
            self.registered_handlers.pop(name)
    
    def list_event_handlers(self):
        return self.registered_handlers.keys()
    
    
        

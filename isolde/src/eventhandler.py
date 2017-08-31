# Copyright 2017 Tristan Croll
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


# vim: set expandtab shiftwidth=4 softtabstop=4:

# Generic wrapper for handling of events. It would also make sense to put
# custom mouse modes here
from collections import defaultdict
class EventHandler():
    def __init__(self, owner):
        self.owner = owner
        self.registered_handlers = {}
        self.handlers_by_trigger_name = defaultdict(list)



    def add_event_handler(self, name, event_type, callback):
        
        handler = self.owner.triggers.add_handler(event_type, callback)
        self.registered_handlers[name] = handler
        self.handlers_by_trigger_name[event_type] = name
        return name


    def remove_event_handler(self, name, error_on_missing = True):
        try:
            handler = self.registered_handlers[name]
            self.owner.triggers.remove_handler(handler)
            self.registered_handlers.pop(name)
        except KeyError as e:
            if error_on_missing:
                raise e

    def remove_all_handlers_for_trigger(self, trigger_name):
        '''
        Remove all handlers that isolde has created for the given trigger.
        '''
        h_list = self.handlers_by_trigger_name[trigger_name]
        for name in h_list:
            self.remove_event_handler(name)
        h_list.clear()

    def remove_all_handlers(self):
        names = list(self.list_event_handlers())
        for name in names:
            handler = self.registered_handlers[name]
            self.owner.triggers.remove_handler(handler)
            self.registered_handlers.pop(name)

    def list_event_handlers(self):
        return self.registered_handlers.keys()

    def list_all_callbacks(self):
        for key, l in self.handlers_by_trigger_name.items():
            print('{}: {}'.format(key, l))

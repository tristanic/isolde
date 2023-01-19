# @Author: Tristan Croll <tic20>
# @Date:   18-Apr-2018
# @Email:  tcroll@altoslabs.com
# @Last modified by:   tic20
# @Last modified time: 26-Apr-2018
# @License: Free for non-commercial use (see license.pdf)
# @Copyright:2016-2019 Tristan Croll



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

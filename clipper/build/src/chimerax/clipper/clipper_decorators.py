from .lib import clipper_python_core as clipper_core

#class MessageStreamSingleton:
    #class __Get_MessageStream:
        #def __init__(self):
            #self.clipper_messages = clipper_core.ClipperMessageStream()
    #instance = None
    #def __init__(self):
        #if not MessageStreamSingleton.instance:
            #MessageStreamSingleton.instance = MessageStreamSingleton.__Get_MessageStream()
    #def __getattr__(self, name):
        #return getattr(self.instance, name)

#_clipper_messages = MessageStreamSingleton().clipper_messages
        



def mappedclass(old_cls):
    '''
     To make things a little more user-friendly than the raw SWIG wrapper,
     all the major Clipper classes are sub-classed here with some extra
     documentation, useful Python methods, @property decorators, etc.
     This requires overcoming one small problem: when Clipper itself creates 
     and returns objects, they're returned as the base class rather than
     the sub-class. So, we have to over-ride the __new__ method for each
     base class to make sure they're instead instantiated as the derived
     class with all the bells and whistles. This function acts as a class
     decorator to automate this. For example, for the class Atom derived
     from clipper_core.Atom, place:
     
     @mappedclass(clipper_core.Atom)
     
     directly before the Atom class declaration.
    '''

    def decorator(cls):
        def __newnew__(this_cls, *args, **kwargs):
            if this_cls == old_cls:
                return object.__new__(cls)
            return object.__new__(this_cls)
        old_cls.__new__ = __newnew__
        
        return cls
    return decorator

def log_clipper(func):
    '''
    Acts as a decorator to direct Clipper messages to the Python console.
    Any messages coming from Clipper are accumulated in _clipper_messages.
    For any core Clipper function which has the potential to generate a
    warning message, simply add the @log_clipper decorator to the Python
    method. Override this function if you want the messages to go somewhere
    else (e.g. to a log file).
    '''
    def func_wrapper(*args, **kwargs):
        _clipper_messages.clear()
        ret = func(*args, **kwargs)
        message_string = _clipper_messages.read_and_clear()
        if message_string:
            print("CLIPPER WARNING:")
            print(message_string)
        return ret
    return func_wrapper


def format_to_string(cls):
    '''
    Class decorator to redirect the Clipper format() function to __str__,
    to provide pretty printing of the object.
    '''
    def format(self):
        return super(self.__class__,self).format()
    def __str__(self):
        return self.format
    setattr(cls, 'format', property(format))
    setattr(cls, '__str__', __str__)
    return cls

def getters_to_properties(*funcs):
    '''
    Class decorator. Add the names of any getter functions with void 
    arguments (e.g. Coord_grid.u()) to convert them to properties. If
    you want the property name to be different from the function name,
    add the desired name and the function name as a tuple 
    (e.g. ('uvw', '_get_uvw'))
    '''
    def property_factory(func):
        def getter(self):
            return getattr(super(self.__class__, self), func)()
        prop = property(getter)
        return prop

    def decorator(cls):
        for func in funcs:
            if type(func) == tuple:
                setattr(cls, func[0], property_factory(func[1]))
            else:
                setattr(cls, func, property_factory(func)) 
        return cls
    return decorator

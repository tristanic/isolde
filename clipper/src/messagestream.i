%inline %{
  namespace clipper
  {
  /*! Clipper messages are directed to an ostream (std::cerr by default).
   * This makes them quite difficult to capture and use sensibly in
   * Python. The solution taken here is to redirect to a stringstream,
   * and apply a simple decorator to key functions in Python to check
   * and report the string buffer as necessary.
   *
   * To make use of this, add code similar to the below to your Python project.
   *
    # You should only ever have a single instance of this object in your project.
    clipper_messages = clipper.ClipperMessageStream()

    def log_clipper(func):
       def func_wrapper(*args, **kwargs):
            clipper_messages.clear()
            func(*args, **kwargs)
            message_string = clipper_messages.read_and_clear()
            if message_string:
                pass # Do whatever you want to do with the messages here
        return func_wrapper
   *
   * For any function with the potential to generate a Clipper warning,
   * add the decorator @log_clipper - i.e.
   *
   * @log_clipper
   * def myfunc(args):
   *    do stuff
   *
   *
   */
  class ClipperMessageStream
  {
    private:
      std::ostringstream stream;
    public:
      ClipperMessageStream()
      {
        this->redirect_clipper_messages();
      }
      ~ClipperMessageStream()
      {
        clipper::Message message;
        message.set_stream(std::cerr);
      }

      std::ostringstream& get_stream ()
      {
        return stream;
      }

      std::ostringstream& redirect_clipper_messages()
      {
        clipper::Message message;
        message.set_stream(stream);
        return stream;
      }

      std::string read_and_clear()
      {
        stream.flush();
        std::string ret = stream.str();
        stream.str("");
        return ret;
      }

      void clear()
      {
        stream.str("");
      }

  };
  };
%}

%pythoncode %{
class MessageStreamSingleton:
  '''
  Creates and maintains a singleton object designed to intercept
  Clipper messages and redirect them from stderr to a string that
  can be used in Python.
  '''
  class __Get_MessageStream:
      def __init__(self):
          self.clipper_messages = ClipperMessageStream()
  instance = None
  def __init__(self):
      if not MessageStreamSingleton.instance:
          MessageStreamSingleton.instance = MessageStreamSingleton.__Get_MessageStream()
  def __getattr__(self, name):
      return getattr(self.instance, name)

_clipper_messages = MessageStreamSingleton().clipper_messages

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


%}


%inline %{
  namespace clipper
  {

  void warn_test()
  {
    Message::message(Message_warn("Test warning"));
  }

  void except_test()
  {
    throw std::length_error("Test exception");
  }
 



  std::string ClipperStringAsString(const clipper::String &a)
  {
    return (std::string)(a);
  }

  }
%}

namespace clipper
{
%pythoncode %{
  warn_test = log_clipper(warn_test)
%}
}

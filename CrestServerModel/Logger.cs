namespace CrestServer
{
    public static class LogManager
    {
        public static Logger GetCurrentClassLogger() => new Logger();
    }
    public class Logger
    {
        public void Error(string wavelengthMustBeF)
        {
            //write to log here
        }
    }
}
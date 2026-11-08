//Config.h
#ifndef CONFIG_H_
#define CONFIG_H_
#include <string>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>


/* \brief Generic configuration Class */
class Config {
    protected:
        typedef std::map<std::string, std::string>::iterator mapi;
        typedef std::map<std::string, std::string>::const_iterator mapci;

        std::string delimiter_;//!< separator between key and value
        std::string comment_;  //!< separator between value and comments
        //!< extracted keys and values
        std::map<std::string,std::string> contents_;

    public:
        Config(std::string filename, std::string delimiter="=",
            std::string comment="#");
        Config();

        // Write or read configuration
        friend std::ostream& operator<<(std::ostream& os, const Config& cf);
        friend std::istream& operator>>(std::istream& is, Config& cf);

        //Search for key and value or optional default value, call as read<T>
        /*template<class T> T Read(const std::string& in_key) const; */
        template<class T> T Read(const std::string& in_key,
            const T& in_value) const;
        template<class T> bool ReadInto(T& out_var,
            const std::string& in_key) const;
        template<class T> bool ReadInto(T& out_var, const std::string& in_key,
            const T& in_value) const;
        bool FileExist(std::string filename);
        void ReadFile(std::string filename, std::string delimiter="=",
            std::string comment="#");

        // Check whether key exists in configuration
        bool KeyExists(const std::string& in_key) const;

        // Modify keys and values
        template<class T> void Add(const std::string& in_key,
            const T& in_value);
        void Remove(const std::string& in_key);

        // Check or change configuration syntax
        std::string GetDelimiter() const;
        std::string GetComment() const;
        std::string SetDelimiter(const std::string& in_s);
        std::string SetComment(const std::string& in_s);

    protected:
        template<class T> static std::string T_as_string(const T& t);
        template<class T> static T string_as_T(const std::string& s);
        static void Trim(std::string& inout_s);

    // Exception types
    public:
        struct File_not_found 
        {
            std::string filename;
            File_not_found(const std::string& filename_=std::string()):
                filename(filename_) {}
        };

        struct Key_not_found 
        {
            // thrown only by T read(key) variant of read()
            std::string key;
            Key_not_found(const std::string& key_=std::string()):
                key(key_) {}
        };
};


Config::Config(std::string filename, std::string delimiter,
    std::string comment): delimiter_(delimiter), comment_(comment)
{
    // Construct a Config, getting keys and values from given file  
    std::ifstream in(filename.c_str());
    if(!in) throw File_not_found(filename);
    in >> (*this);
};

Config::Config(): delimiter_(std::string(1,'=')), comment_(std::string(1,'#'))
{
    // Construct a Config without a file; empty
};

std::ostream& operator<<(std::ostream& os, const Config& cf)
{
    // Save a Config to os
    for(Config::mapci p = cf.contents_.begin(); p != cf.contents_.end(); ++p)
    {
        os << p->first << " " << cf.delimiter_ << " ";
        os << p->second << std::endl;
    }
    return os;
};

std::istream& operator>>(std::istream& is, Config& cf)
{
    // Load a Config from is
    // Read in keys and values, keeping internal whitespace
    typedef std::string::size_type pos;
    const std::string& delim  = cf.delimiter_;// separator
    const std::string& comm   = cf.comment_;  // comment
    const pos skip = delim.length();      // length of separator

    // might need to read ahead to see where value ends
    std::string nextline = "";

    while(is || nextline.length() > 0)
    {
        // Read an entire line at a time
        std::string line;
        if(nextline.length() > 0)
        {
            line = nextline;// we read ahead; use it now  
            nextline = "";
        }
        else
        {
            std::getline(is, line);
        }

        // Ignore comments
        line = line.substr(0, line.find(comm));

        // Parse the line if it contains a delimiter
        pos delimPos = line.find(delim);
        if(delimPos < std::string::npos)
        {
            // Extract the key  
            std::string key = line.substr(0, delimPos);
            line.replace(0, delimPos+skip, "");

            // See if value continues on the next line
            // Stop at blank line, next line with a key, end of stream,
            // or end of file sentry  
            bool terminate = false;
            while(!terminate && is)
            {
                std::getline(is, nextline);
                terminate = true;

                std::string nlcopy = nextline;
                Config::Trim(nlcopy);
                if(nlcopy == "") continue;

                nextline = nextline.substr(0, nextline.find(comm));
                if(nextline.find(delim) != std::string::npos)
                    continue;

                nlcopy = nextline;
                Config::Trim(nlcopy);
                if(nlcopy != "") line += "\n";
                line += nextline;
                terminate = false;
            }

            // Store key and value
            Config::Trim(key);
            Config::Trim(line);
            cf.contents_[key] = line;// overwrites if key is repeated  
        }
    }

    return is;
};

/*
template<class T> T Config::Read(const std::string& key) const
{
    // Read the value corresponding to key
    mapci p = contents_.find(key);
    if(p == contents_.end()) throw Key_not_found(key);
    return string_as_T<T>(p->second);
};
*/

template<class T> T Config::Read(const std::string& key, const T& value) const
{
    // Return the value corresponding to key or given default value
    // if key is not found
    mapci p = contents_.find(key);
    // if(p == contents_.end()) return value;
    if (p == contents_.end()) throw Key_not_found(key);
    return string_as_T<T>(p->second);
};

template<class T> bool Config::ReadInto(T& var, const std::string& key) const
{
    // Get the value corresponding to key and store in var
    // Return true if key is found
    // Otherwise leave var untouched
    mapci p = contents_.find(key);
    bool found = (p != contents_.end());
    if(found) var = string_as_T<T>(p->second);
    return found;
};

template<class T> bool Config::ReadInto(T& var, const std::string& key,
    const T& value) const
{
    // Get the value corresponding to key and store in var
    // Return true if key is found
    // Otherwise set var to given default
    mapci p = contents_.find(key);
    bool found = (p != contents_.end());
    if(found)
        var = string_as_T<T>(p->second);
    else
        var = value;
    return found;
};

bool Config::FileExist(std::string filename)
{
    bool exist= false;
    std::ifstream in(filename.c_str());
    if(in) exist = true;
    return exist;
};

void Config::ReadFile(std::string filename, std::string delimiter,
    std::string comment)
{
    delimiter_ = delimiter;
    comment_ = comment;
    std::ifstream in(filename.c_str());
    if(!in) throw File_not_found(filename);
    in >> (*this);
};

bool Config::KeyExists(const std::string& key) const
{
    // Indicate whether key is found
    mapci p = contents_.find(key);
    return (p != contents_.end());
};

template<class T> void Config::Add(const std::string& in_key, const T& value)
{
    // Add a key with given value
    std::string v = T_as_string(value);
    std::string key=in_key;
    Trim(key);
    Trim(v);
    contents_[key] = v;
    return;
};

void Config::Remove(const std::string& key)
{
    // Remove key and its value
    contents_.erase(contents_.find(key));
    return;
};

// Check or change configuration syntax
std::string Config::GetDelimiter() const
{
    return delimiter_;
};

std::string Config::GetComment() const 
{
    return comment_;
};

std::string Config::SetDelimiter(const std::string& in_s)
{
    std::string old = delimiter_;delimiter_ = in_s;return old;
};

std::string Config::SetComment(const std::string& in_s)
{
    std::string old = comment_;comment_ =  in_s;return old;
};

/* static */
template<class T> std::string Config::T_as_string(const T& t)
{
    // Convert from a T to a string
    // Type T must support << operator
    std::ostringstream ost;
    ost << t;
    return ost.str();
};

/* static */
template<class T> T Config::string_as_T(const std::string& s)
{
    // Convert from a string to a T
    // Type T must support >> operator
    T t;
    std::istringstream ist(s);
    ist >> t;
    return t;
};

/* static */
template<>
inline std::string Config::string_as_T<std::string>(const std::string& s)
{
    // Convert from a string to a string
    // In other words, do nothing
    return s;
};

/* static */
template<>
inline bool Config::string_as_T<bool>(const std::string& s)  
{
    // Convert from a string to a bool
    // Interpret "false", "F", "no", "n", "0" as false  
    // Interpret "true", "T", "yes", "y", "1", "-1", or anything else as true  
    bool b = true;
    std::string sup = s;
    for(std::string::iterator p = sup.begin(); p != sup.end(); ++p)  
        *p = toupper(*p);// make string all caps  
    if(sup==std::string("FALSE") || sup==std::string("F") ||  
        sup==std::string("NO") || sup==std::string("N") ||  
        sup==std::string("0") || sup==std::string("NONE"))  
        b = false;
    return b;
};

/* static */
void Config::Trim(std::string& inout_s)
{
    // Remove leading and trailing whitespace
    static const char whitespace[] = " \n\t\v\r\f";
    inout_s.erase(0, inout_s.find_first_not_of(whitespace));
    inout_s.erase(inout_s.find_last_not_of(whitespace) + 1U);
};

#endif
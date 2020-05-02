#ifndef JSONSTREAM_H
#define JSONSTREAM_H

#include <fstream>
#include <iostream>
#include <regex>

#include "DisambiguateSTLContainer.h"



//! Represents a \c string in json-format, i.e. with leading and trailing '"'
class JsonString {
public:
    //! c'tor (1) to create a json string from a std::string
    JsonString(std::string const & string) : value_{} {
        createValidString(string);
    }
    //! c'tor (2) to create a json string from any other valid type
    template<typename ValueType>
    JsonString(ValueType const & value) : value_{} {
        createValidString(std::to_string(value));
    }

    //! Getter for member \c value_
    auto const & value() const { return value_; }

private:
    void createValidString(std::string string) {
        if (string.size() == 0) {
            value_ = "\"\"";    // ""
        } else if (string.size() == 1 && string == "\"") {
            value_ = "\"\\\"\"";    // "\""
        } else {
            // remove trailing ", re-add later
            if (string.at(0) == '"' && string.at(string.size()-1) == '"') {
                string.replace(0,1,"");
                string.replace(string.size()-1,1,"");
            }
            size_t nextSearchPos = 0;
            auto pos = string.find('"', nextSearchPos);
            while (pos != std::string::npos) {
                nextSearchPos = pos;
                if (pos == 0 || (pos > 0 && string.at(pos-1) != '\\')) { // unescaped " at pos
                    string.replace(pos, 1, "\\\""); // Replace " by \"
                    nextSearchPos += 2;
                } else {
                    ++nextSearchPos;
                }
                pos = string.find('"', nextSearchPos);
            }
        }
        value_ = std::string("\"") + string + std::string("\"");    // add " at beginning and end
    }

    //! Stores json string
    std::string value_;
};



//! Represents a \c value in json-format
/*! Can create (nested) json representations of
 * string, number, object, array, true and false.
 * 'null' is not implemented (yet) */
class JsonValue {
public:
    // tag dispatch for explicitly turn a string into a JsonValue without modification
    static struct StringIsValidJson{} stringIsValidJson;

    //! default c'tor, initializes empty \c value_ string
    JsonValue() : value_{} {}
    //! c'tor (1) to create a JSON string from a std::string
    JsonValue(std::string const & string) : value_{JsonString(string).value()} {}
    //! c'tor (2) to turn a std::string into a JsonValue without modification
    JsonValue(std::string const & string, StringIsValidJson) : value_{string} {}
    //! c'tor (3) to create a JSON bool string from a bool
    JsonValue(bool b) : value_{b ? "true" : "false"} {}
    //! c'tor (4) to create a JSON string from any other valid type
    template<typename ValueType>
    JsonValue(ValueType const & value) : value_{} {
        createValue<ValueType>(value, is_container(value));
    }

    //! Compares the string representations of the values by string operator ==
    bool operator==(JsonValue const & rhs) const {
        return value_ == rhs.value_;
    }
    //! Compares the string representations of the values by string operator <
    bool operator<(JsonValue const & rhs) const {
        return value_ < rhs.value_;
    }
    friend std::ostream & operator<<(std::ostream & out, JsonValue const & v) {
        out << v.value_;
        return out;
    }
    //! Getter for member \c value_
    auto const & value() const { return value_; }

private:
    //! Create a JSON string from a sequential container like std::set or std::vector
    template<typename Container>    // sequential container
    void createContainerRepresentation(Container const & container, std::false_type) {
        value_ = std::string("[");
        if (container.size()) {
           for (auto&& elem : container) {
               auto val = JsonValue(elem);
               value_ += val.value();
               value_ += ",";
           }
           *(value_.rbegin()) = ']'; // replace last ',' by ']'
        } else {
            value_ += "]";
        }
    }
    //! Create a JSON string from an associative container like std::map
    template<typename Container>    // associative container
    void createContainerRepresentation(Container const & container, std::true_type) {
        value_ = std::string("{");
        if (container.size()) {
            for (auto&& elem : container) {
                auto key = elem.first;
                auto jkey = JsonString(key).value();
                auto val = elem.second;
                auto jval = JsonValue(val).value();
                value_ += jkey;
                value_ += ":";
                value_ += jval;
                value_ += ",";
            }
            *(value_.rbegin()) = '}'; // replace last ',' by '}'
        } else {
            value_ += "}";
        }
    }
    //! Create a JSON string from a numeric value
    template<typename Number>
    void createNumericRepresentation(Number num) {
        value_ = std::to_string(num);
    }
    //! Create JSON string from any container
    template<typename ValueType>    // container
    void createValue(ValueType const & value, std::true_type) {
        createContainerRepresentation<ValueType>(value, is_map_container(value));
    }
    //! Create JSON string from a non-container value
    template<typename ValueType>    // number
    void createValue(ValueType const & value, std::false_type) {
        createNumericRepresentation<ValueType>(value);
    }

    //! Stores json-representation of input value
    std::string value_;
};



//! Write a large array or object in json representation to a file
class JsonStream {
public:
    JsonStream(std::basic_ostream<char> & outstream)
        : outstream_{outstream} {}
    auto& ostream() const { return outstream_; }
protected:
    std::basic_ostream<char> & outstream_;
};



//! Write elements of an array to a json representation in a file
class JsonStreamArray : public JsonStream {
public:
    //! Begin a new json array, write immediately into \c outstream
    JsonStreamArray(std::basic_ostream<char> & outstream, bool newline = false)
        : JsonStream(outstream),
          open_{true}, newline_{newline}, setKomma_{false} {
        this->outstream_ << "[";
    }
    //! If json array was not closed, do this on destruction
    ~JsonStreamArray() {
        if (open_) {
            close();
        }
    }
    //! Append new JsonValue to the array
    void addValue(JsonValue const & value) {
        if (!open_) { throw std::runtime_error("[ERROR] -- JsonStreamArray -- Tried to write to a closed JsonStream"); }
        if (setKomma_) {
            this->outstream_ << ",";
            if (newline_) { this->outstream_ << std::endl; }
        }
        this->outstream_ << value.value();
        setKomma_ = true;
    }
    //! Append a new element to the array from any valid type
    template<typename ValueType>
    void addValue(ValueType const & value) {
        addValue(JsonValue(value));
    }
    //! Append a new element to the array using stream operator
    template<typename ValueType>
    JsonStreamArray & operator<<(ValueType const & value) {
        addValue(value);
        return *this;
    }
    //! Close json array
    void close() {
        this->outstream_ << "]";
        open_ = false;
        setKomma_ = false;
    }
private:
    //! Flag if json array is open (i.e. no closing "]" was written)
    bool open_;
    //! Flag if insert newline after each element
    bool newline_;
    //! Flag if need to set a komma before a new element
    bool setKomma_;
};



//! Write elements of a dict to a json representation in a file
class JsonStreamDict : public JsonStream {
public:
    //! Begin a new json dict, write immediately into \c outstream
    JsonStreamDict(std::basic_ostream<char> & outstream, bool newline = false)
        : JsonStream(outstream),
          open_{true}, newline_{newline}, setKomma_{false} {
        this->outstream_ << "{";
    }
    //! If json dict was not closed, do this on destruction
    ~JsonStreamDict() {
        if (open_) {
            close();
        }
    }
    //! Append a new key: JsonValue pair to the dict
    void addValue(std::string const & key, JsonValue const & value) {
        if (!open_) { throw std::runtime_error("[ERROR] -- JsonStreamArray -- Tried to write to a closed JsonStream"); }
        if (setKomma_) {
            this->outstream_ << ",";
            if (newline_) { this->outstream_ << std::endl; }
        }
        this->outstream_ << JsonValue(key).value() << ":" << value.value();
        setKomma_ = true;
    }
    //! Append a new key:value pair to the dict
    template<typename ValueType>
    void addValue(std::string const & key, ValueType const & value) {
        addValue(key, JsonValue(value));
    }
    //! Close json dict
    void close() {
        this->outstream_ << "}";
        open_ = false;
        setKomma_ = false;
    }
private:
    //! Flag if json dict is open (i.e. no closing "}" was written)
    bool open_;
    //! Flag if insert newline after each key:value pair
    bool newline_;
    //! Flag if need to set a komma before a new key:value pair
    bool setKomma_;
};

#endif // JSONSTREAM_H

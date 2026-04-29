#include <string>
#include <regex>
#include <limits.h>

std::string decodeBits(const std::string& bits) {
    auto start = bits.find_first_of('1');
    auto end = bits.find_last_of('1');
    std::string result = bits.substr(start, end-start+1);
    std::smatch sm;
    std::string str = result;
    int unit = INT_MAX;
    while (std::regex_search(str, sm, std::regex("0+|1+")) ) {
        if( sm[0].length() < unit )
            unit = sm[0].length();
        str = sm.suffix();
    }
    while( result.length() % unit != 0 ) unit /= 3;
    std::regex e("0{"+std::to_string(unit)+"}");
    result = std::regex_replace(result, e, "0");
    e = std::regex("1{"+std::to_string(unit)+"}");
    result = std::regex_replace(result, e, "1");
    result = std::regex_replace(result, std::regex("0000000"), "   ");
    result = std::regex_replace(result, std::regex("111"), "-");
    result = std::regex_replace(result, std::regex("000"), " ");
    result = std::regex_replace(result, std::regex("1"), ".");
    result = std::regex_replace(result, std::regex("0"), "");
    return result;
}

std::string decodeMorse(const std::string& morse) {
    std::string decoded = "";
    std::string word = "";
    int nbSpace = 0;
    bool isStarted = false;
    for( auto p : morse ) {
      if ( p != ' ' ) {
        if ( nbSpace == 3 )
          decoded += " ";
        word += p;
        nbSpace = 0;
        isStarted = true;
      }
      else if ( isStarted ) {
        nbSpace++;
        if ( nbSpace == 1 ) {
          decoded += MORSE_CODE[word];
          word = "";
        }
      }       
    }
    decoded += MORSE_CODE[word];
    return decoded;
}

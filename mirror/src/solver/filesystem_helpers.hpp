#pragma once

#include <string>

namespace openvfife {

// Legacy filesystem helpers remain private until #46 isolates result I/O.
bool isFileExists(const std::string& fname);
bool isDirExists(const std::string& path);
bool makedir(const std::string& path);

} // namespace openvfife

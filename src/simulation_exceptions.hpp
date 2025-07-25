#pragma once

#include <exception>
#include <string>

class SimulationException : public std::exception {
protected:
    std::string message;
    std::string error_type;

public:
    SimulationException(const std::string& msg, const std::string& type) 
        : message(msg), error_type(type) {}
    
    const char* what() const noexcept override {
        return message.c_str();
    }
    
    const std::string& get_error_type() const noexcept {
        return error_type;
    }
};

class ConvergenceException : public SimulationException {
public:
    ConvergenceException(const std::string& msg, const std::string& type = "CONVERGENCE_FAILURE") 
        : SimulationException(msg, type) {}
};

class SetupException : public SimulationException {
public:
    SetupException(const std::string& msg, const std::string& type = "SETUP_FAILURE") 
        : SimulationException(msg, type) {}
};

class MatrixException : public SimulationException {
public:
    MatrixException(const std::string& msg, const std::string& type = "MATRIX_FAILURE") 
        : SimulationException(msg, type) {}
};

class DeviceException : public SimulationException {
public:
    DeviceException(const std::string& msg, const std::string& type = "DEVICE_FAILURE") 
        : SimulationException(msg, type) {}
};

class ParsingException : public SimulationException {
public:
    ParsingException(const std::string& msg, const std::string& type = "PARSING_FAILURE") 
        : SimulationException(msg, type) {}
};
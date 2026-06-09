#ifndef RESULTS_STORAGE_H
#define RESULTS_STORAGE_H

#include <Eigen/Dense>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <H5Cpp.h>

using json = nlohmann::json;

/**
 * @brief Helper class to store simulation results in HDF5 format.
 *
 * Provides functions to write Eigen vectors and matrices to HDF5 files
 * in a format compatible with the Python ResultsStorage.load() function.
 */
class ResultsStorage
{
public:
    /**
     * @brief Create a new HDF5 file for storing results
     * @param filename Path to the output HDF5 file
     */
    ResultsStorage(const std::string &filename, bool overwrite = false);

    ~ResultsStorage();

    /**
     * @brief Open an existing HDF5 file for reading
     * @param filename Path to the HDF5 file to open
     * @return ResultsStorage instance with file opened in read mode
     */
    static ResultsStorage openForReading(const std::string &filename);

    // Metadata and configuration methods
    void setSolverSettings(const json &settings);
    void setModelSettings(const json &settings);
    void setStorageConfig(bool energy, bool power, bool drift,
                          bool sav, bool solver_setting, bool model_setting);
    void setSuccess(bool success);

    // Methods to write 1D arrays (vectors)
    void writeVector(const std::string &name, const Eigen::VectorXd &vec);
    void writeVector(const std::string &name, const Eigen::VectorXf &vec);
    void writeVector(const std::string &name, const std::vector<double> &vec);
    void writeVector(const std::string &name, const std::vector<float> &vec);

    // Methods to write 2D arrays (matrices)
    void writeMatrix(const std::string &name, const Eigen::MatrixXd &mat);
    void writeMatrix(const std::string &name, const Eigen::MatrixXf &mat);

    // Methods to write scalar values as attributes
    void writeAttribute(const std::string &name, double value);
    void writeAttribute(const std::string &name, float value);
    void writeAttribute(const std::string &name, int value);
    void writeAttribute(const std::string &name, bool value);
    void writeAttribute(const std::string &name, const std::string &value);
    void writeAttribute(const std::string &name, const json &value);

    // Convenience method to write indices
    void writeIndices(const std::string &name, const Eigen::VectorXi &indices);
    void writeIndices(const std::string &name, const std::vector<int> &indices);

    // Methods to read attributes
    bool readAttribute(const std::string &name, double &value);
    bool readAttribute(const std::string &name, float &value);
    bool readAttribute(const std::string &name, int &value);
    bool readAttribute(const std::string &name, bool &value);
    bool readAttribute(const std::string &name, std::string &value);
    bool readAttribute(const std::string &name, json &value);

    // Methods to read 1D arrays (vectors)
    bool readVector(const std::string &name, Eigen::VectorXd &vec);
    bool readVector(const std::string &name, Eigen::VectorXf &vec);
    bool readVector(const std::string &name, std::vector<double> &vec);
    bool readVector(const std::string &name, std::vector<float> &vec);

    // Methods to read 2D arrays (matrices)
    bool readMatrix(const std::string &name, Eigen::MatrixXd &mat);
    bool readMatrix(const std::string &name, Eigen::MatrixXf &mat);

    // Methods to read indices
    bool readIndices(const std::string &name, Eigen::VectorXi &indices);
    bool readIndices(const std::string &name, std::vector<int> &indices);

    // Check if a dataset or attribute exists
    bool exists(const std::string &name);

    // Close and finalize the file
    void close();

private:
    std::string filename_;
    H5::H5File *file_;
    bool is_open_;

    json storage_config_;
};

#endif // RESULTS_STORAGE_H

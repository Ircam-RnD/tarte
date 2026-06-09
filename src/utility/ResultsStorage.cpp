#include "ResultsStorage.h"
#include <iostream>
#include <filesystem>

ResultsStorage::ResultsStorage(const std::string &filename, bool overwrite)
    : filename_(filename), file_(nullptr), is_open_(false)
{
    try
    {
        // Create parent directory if it doesn't exist
        std::filesystem::path filepath(filename);
        if (filepath.has_parent_path())
        {
            std::filesystem::create_directories(filepath.parent_path());
        }

        // Create or open the HDF5 file (without overwriting if exists)
        H5::FileAccPropList plist;
        if (overwrite)
        {
            file_ = new H5::H5File(filename, H5F_ACC_TRUNC, H5::FileCreatPropList::DEFAULT, plist);
            is_open_ = true;
        }
        else
        {
            file_ = new H5::H5File(filename, H5F_ACC_CREAT | H5F_ACC_RDWR, H5::FileCreatPropList::DEFAULT, plist);
            is_open_ = true;
        }
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error creating HDF5 file: " << e.getCDetailMsg() << std::endl;
        is_open_ = false;
    }
}

ResultsStorage::~ResultsStorage()
{
    close();
}

void ResultsStorage::setStorageConfig(bool energy, bool power, bool drift,
                                      bool sav, bool solver_setting, bool model_setting)
{
    storage_config_["Energy"] = energy;
    storage_config_["Power"] = power;
    storage_config_["Drift"] = drift;
    storage_config_["SAV"] = sav;
    storage_config_["SolverSetting"] = solver_setting;
    storage_config_["ModelSetting"] = model_setting;
}

void ResultsStorage::setSolverSettings(const json &settings)
{
    if (is_open_ && file_)
    {
        try
        {
            H5::StrType str_type(H5::PredType::C_S1, settings.dump().size());
            H5::Attribute attr = file_->createAttribute("solver_dict", str_type, H5S_SCALAR);
            attr.write(str_type, settings.dump());
        }
        catch (const H5::Exception &e)
        {
            std::cerr << "Error writing solver_dict: " << e.getCDetailMsg() << std::endl;
        }
    }
}

void ResultsStorage::setModelSettings(const json &settings)
{
    if (is_open_ && file_)
    {
        try
        {
            H5::StrType str_type(H5::PredType::C_S1, settings.dump().size());
            H5::Attribute attr = file_->createAttribute("model_dict", str_type, H5S_SCALAR);
            attr.write(str_type, settings.dump());
        }
        catch (const H5::Exception &e)
        {
            std::cerr << "Error writing model_dict: " << e.getCDetailMsg() << std::endl;
        }
    }
}

void ResultsStorage::setSuccess(bool success)
{
    if (is_open_ && file_)
    {
        try
        {
            std::string success_str = json(success).dump();
            H5::StrType str_type(H5::PredType::C_S1, success_str.size());
            H5::Attribute attr = file_->createAttribute("success", str_type, H5S_SCALAR);
            attr.write(str_type, success_str);
        }
        catch (const H5::Exception &e)
        {
            std::cerr << "Error writing success: " << e.getCDetailMsg() << std::endl;
        }
    }
}

void ResultsStorage::writeVector(const std::string &name, const Eigen::VectorXd &vec)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        hsize_t dims[1] = {(hsize_t)vec.size()};
        H5::DataSpace dataspace(1, dims);
        H5::DataSet dataset = file_->createDataSet(name, H5::PredType::IEEE_F64LE, dataspace);
        dataset.write(vec.data(), H5::PredType::NATIVE_DOUBLE);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing vector '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeVector(const std::string &name, const Eigen::VectorXf &vec)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        hsize_t dims[1] = {(hsize_t)vec.size()};
        H5::DataSpace dataspace(1, dims);
        H5::DataSet dataset = file_->createDataSet(name, H5::PredType::IEEE_F32LE, dataspace);
        dataset.write(vec.data(), H5::PredType::NATIVE_FLOAT);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing vector '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeVector(const std::string &name, const std::vector<double> &vec)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        hsize_t dims[1] = {(hsize_t)vec.size()};
        H5::DataSpace dataspace(1, dims);
        H5::DataSet dataset = file_->createDataSet(name, H5::PredType::IEEE_F64LE, dataspace);
        dataset.write(vec.data(), H5::PredType::NATIVE_DOUBLE);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing vector '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeVector(const std::string &name, const std::vector<float> &vec)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        hsize_t dims[1] = {(hsize_t)vec.size()};
        H5::DataSpace dataspace(1, dims);
        H5::DataSet dataset = file_->createDataSet(name, H5::PredType::IEEE_F32LE, dataspace);
        dataset.write(vec.data(), H5::PredType::NATIVE_FLOAT);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing vector '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeMatrix(const std::string &name, const Eigen::MatrixXd &mat)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        hsize_t dims[2] = {(hsize_t)mat.rows(), (hsize_t)mat.cols()};
        H5::DataSpace dataspace(2, dims);
        H5::DataSet dataset = file_->createDataSet(name, H5::PredType::IEEE_F64LE, dataspace);
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat_rowmajor = mat;
        dataset.write(mat_rowmajor.data(), H5::PredType::NATIVE_DOUBLE);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing matrix '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeMatrix(const std::string &name, const Eigen::MatrixXf &mat)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        hsize_t dims[2] = {(hsize_t)mat.rows(), (hsize_t)mat.cols()};
        H5::DataSpace dataspace(2, dims);
        H5::DataSet dataset = file_->createDataSet(name, H5::PredType::IEEE_F32LE, dataspace);
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> mat_rowmajor = mat;
        dataset.write(mat_rowmajor.data(), H5::PredType::NATIVE_FLOAT);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing matrix '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeIndices(const std::string &name, const Eigen::VectorXi &indices)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        hsize_t dims[1] = {(hsize_t)indices.size()};
        H5::DataSpace dataspace(1, dims);
        H5::DataSet dataset = file_->createDataSet(name, H5::PredType::STD_I32LE, dataspace);
        dataset.write(indices.data(), H5::PredType::NATIVE_INT);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing indices '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeIndices(const std::string &name, const std::vector<int> &indices)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        hsize_t dims[1] = {(hsize_t)indices.size()};
        H5::DataSpace dataspace(1, dims);
        H5::DataSet dataset = file_->createDataSet(name, H5::PredType::STD_I32LE, dataspace);
        dataset.write(indices.data(), H5::PredType::NATIVE_INT);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing indices '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeAttribute(const std::string &name, double value)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        H5::Attribute attr = file_->createAttribute(name, H5::PredType::IEEE_F64LE, H5S_SCALAR);
        attr.write(H5::PredType::NATIVE_DOUBLE, &value);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing attribute '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeAttribute(const std::string &name, float value)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        H5::Attribute attr = file_->createAttribute(name, H5::PredType::IEEE_F32LE, H5S_SCALAR);
        attr.write(H5::PredType::NATIVE_FLOAT, &value);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing attribute '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeAttribute(const std::string &name, int value)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        H5::Attribute attr = file_->createAttribute(name, H5::PredType::STD_I32LE, H5S_SCALAR);
        attr.write(H5::PredType::NATIVE_INT, &value);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing attribute '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeAttribute(const std::string &name, bool value)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        std::string bool_str = json(value).dump();
        H5::StrType str_type(H5::PredType::C_S1, bool_str.size());
        H5::Attribute attr = file_->createAttribute(name, str_type, H5S_SCALAR);
        attr.write(str_type, bool_str);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing attribute '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeAttribute(const std::string &name, const std::string &value)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        H5::StrType str_type(H5::PredType::C_S1, value.size());
        H5::Attribute attr = file_->createAttribute(name, str_type, H5S_SCALAR);
        attr.write(str_type, value);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing attribute '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::writeAttribute(const std::string &name, const json &value)
{
    if (!is_open_ || !file_)
        return;

    try
    {
        std::string json_str = value.dump();
        H5::StrType str_type(H5::PredType::C_S1, json_str.size());
        H5::Attribute attr = file_->createAttribute(name, str_type, H5S_SCALAR);
        attr.write(str_type, json_str);
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error writing attribute '" << name << "': " << e.getCDetailMsg() << std::endl;
    }
}

void ResultsStorage::close()
{
    if (is_open_ && file_)
    {
        try
        {
            // Write storage configuration if set
            if (!storage_config_.empty())
            {
                std::string config_str = storage_config_.dump();
                H5::StrType str_type(H5::PredType::C_S1, config_str.size());
                H5::Attribute attr = file_->createAttribute("storage_config", str_type, H5S_SCALAR);
                attr.write(str_type, config_str);
            }

            file_->close();
            delete file_;
            file_ = nullptr;
            is_open_ = false;
        }
        catch (const H5::Exception &e)
        {
            std::cerr << "Error closing HDF5 file: " << e.getCDetailMsg() << std::endl;
        }
    }
}

// Static factory method to open existing files for reading
ResultsStorage ResultsStorage::openForReading(const std::string &filename)
{
    ResultsStorage storage(filename);
    if (storage.file_)
    {
        try
        {
            storage.file_->close();
            delete storage.file_;
            storage.file_ = nullptr;
        }
        catch (...)
        {
        }
    }

    try
    {
        storage.file_ = new H5::H5File(filename, H5F_ACC_RDONLY);
        storage.is_open_ = true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error opening HDF5 file for reading: " << e.getCDetailMsg() << std::endl;
        storage.is_open_ = false;
    }
    return storage;
}

// Check if a dataset or attribute exists
bool ResultsStorage::exists(const std::string &name)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        return file_->nameExists(name);
    }
    catch (...)
    {
        return false;
    }
}

// Read attribute methods
bool ResultsStorage::readAttribute(const std::string &name, double &value)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::Attribute attr = file_->openAttribute(name);
        attr.read(H5::PredType::NATIVE_DOUBLE, &value);
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading attribute '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

bool ResultsStorage::readAttribute(const std::string &name, float &value)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::Attribute attr = file_->openAttribute(name);
        attr.read(H5::PredType::NATIVE_FLOAT, &value);
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading attribute '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

bool ResultsStorage::readAttribute(const std::string &name, int &value)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::Attribute attr = file_->openAttribute(name);
        attr.read(H5::PredType::NATIVE_INT, &value);
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading attribute '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

bool ResultsStorage::readAttribute(const std::string &name, bool &value)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::Attribute attr = file_->openAttribute(name);
        H5::StrType str_type = attr.getStrType();
        std::string str_val;
        char *buf = new char[attr.getStorageSize() + 1];
        attr.read(str_type, buf);
        str_val = std::string(buf);
        delete[] buf;

        json j = json::parse(str_val);
        value = j.get<bool>();
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading attribute '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

bool ResultsStorage::readAttribute(const std::string &name, std::string &value)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::Attribute attr = file_->openAttribute(name);
        H5::StrType str_type = attr.getStrType();
        std::string result(attr.getStorageSize(), '\0');
        attr.read(str_type, &result[0]);
        value = result;
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading attribute '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

bool ResultsStorage::readAttribute(const std::string &name, json &value)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::Attribute attr = file_->openAttribute(name);
        H5::StrType str_type = attr.getStrType();
        std::string str_val(attr.getStorageSize(), '\0');
        attr.read(str_type, &str_val[0]);
        value = json::parse(str_val);
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading attribute '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

// Read vector methods
bool ResultsStorage::readVector(const std::string &name, Eigen::VectorXd &vec)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::DataSet dataset = file_->openDataSet(name);
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dims[1];
        dataspace.getSimpleExtentDims(dims, nullptr);

        vec.resize(dims[0]);
        dataset.read(vec.data(), H5::PredType::NATIVE_DOUBLE);
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading vector '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

bool ResultsStorage::readVector(const std::string &name, Eigen::VectorXf &vec)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::DataSet dataset = file_->openDataSet(name);
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dims[1];
        dataspace.getSimpleExtentDims(dims, nullptr);

        vec.resize(dims[0]);
        dataset.read(vec.data(), H5::PredType::NATIVE_FLOAT);
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading vector '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

bool ResultsStorage::readVector(const std::string &name, std::vector<double> &vec)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::DataSet dataset = file_->openDataSet(name);
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dims[1];
        dataspace.getSimpleExtentDims(dims, nullptr);

        vec.resize(dims[0]);
        dataset.read(vec.data(), H5::PredType::NATIVE_DOUBLE);
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading vector '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

bool ResultsStorage::readVector(const std::string &name, std::vector<float> &vec)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::DataSet dataset = file_->openDataSet(name);
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dims[1];
        dataspace.getSimpleExtentDims(dims, nullptr);

        vec.resize(dims[0]);
        dataset.read(vec.data(), H5::PredType::NATIVE_FLOAT);
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading vector '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

// Read matrix methods
bool ResultsStorage::readMatrix(const std::string &name, Eigen::MatrixXd &mat)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::DataSet dataset = file_->openDataSet(name);
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dims[2];
        dataspace.getSimpleExtentDims(dims, nullptr);

        mat.resize(dims[0], dims[1]);
        dataset.read(mat.data(), H5::PredType::NATIVE_DOUBLE);
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading matrix '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

bool ResultsStorage::readMatrix(const std::string &name, Eigen::MatrixXf &mat)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::DataSet dataset = file_->openDataSet(name);
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dims[2];
        dataspace.getSimpleExtentDims(dims, nullptr);

        mat.resize(dims[0], dims[1]);
        dataset.read(mat.data(), H5::PredType::NATIVE_FLOAT);
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading matrix '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

// Read indices methods
bool ResultsStorage::readIndices(const std::string &name, Eigen::VectorXi &indices)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::DataSet dataset = file_->openDataSet(name);
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dims[1];
        dataspace.getSimpleExtentDims(dims, nullptr);

        indices.resize(dims[0]);
        dataset.read(indices.data(), H5::PredType::NATIVE_INT);
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading indices '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

bool ResultsStorage::readIndices(const std::string &name, std::vector<int> &indices)
{
    if (!is_open_ || !file_)
        return false;

    try
    {
        H5::DataSet dataset = file_->openDataSet(name);
        H5::DataSpace dataspace = dataset.getSpace();
        hsize_t dims[1];
        dataspace.getSimpleExtentDims(dims, nullptr);

        indices.resize(dims[0]);
        dataset.read(indices.data(), H5::PredType::NATIVE_INT);
        return true;
    }
    catch (const H5::Exception &e)
    {
        std::cerr << "Error reading indices '" << name << "': " << e.getCDetailMsg() << std::endl;
        return false;
    }
}

///
/// \class femtofitter/src/DataIterator3D.hpp
///

#include "Data3D.hpp"


/// \class DataIterator3D
/// \brief Abstract the behavior of looping over data
///
struct DataIterator3D {

  std::vector<Data3D::Datum>::const_iterator
    it,
    end;

  DataIterator3D(const Data3D &data)
    : it(data.data.cbegin())
    , end(data.data.cend())
    {}

  /// Returns true if datum is filled
  bool Next(Data3D::Datum &dat)
    {
      if (it == end || ++it == end) {
        return false;
      }

      dat = *it;
      return true;
    }

  const Data3D::Datum& Get() const
    {
      return *it;
    }
};

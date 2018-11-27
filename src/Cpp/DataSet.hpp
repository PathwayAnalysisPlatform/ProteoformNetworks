#include <string>
#include <set>
#include <EntitySet.hpp>

class DataSet {
  public:
   DataSet();
   ~DataSet();
   void loadData();
   private:
    std::string name;
    std::set<EntitySet> sets;

};
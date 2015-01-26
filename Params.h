#include <iostream>
#include <string.h> // memset
#include <stdio.h>
#include <iomanip>      // std::setprecision
#include <sqlite3pp.h>

using namespace std;

class Params{
public:
    int const N ;
    int const T ;
    const size_t GRID_X_POINTS ;

    size_t dimensionality = 1;
    float D ;
    float dt;
    float k;
    int const RANDOMNESS_BITS_U = 10;
    float VARIANCE;
 
//    const char * RUN_NAME;
//    const char * OUT_FILE ;
    string OUT_FILE_STR ;
    string RUN_NAME_STR ;

    Params(int n, int t, int gp, float d, float k_, float dt_): N(n), T(t), GRID_X_POINTS(gp), D(d), k(k_), dt(dt_)
    {
        VARIANCE = 2.0 * (float) dimensionality * D;
    }

    Params(int n, int t, int gp, float d, float k_, float dt_, const char * rn): \
    N(n), T(t), GRID_X_POINTS(gp), D(d), k(k_), dt(dt_), RUN_NAME_STR(rn)
    {
        VARIANCE = 2.0 * (float) dimensionality * D;
    }

    Params(int n, int t, int gp, float d, float k_, float dt_, const char * rn, const char * out): \
    N(n), T(t), GRID_X_POINTS(gp), D(d), k(k_), dt(dt_), RUN_NAME_STR(rn), OUT_FILE_STR(out)
    {
//        const char * OUT_FILE = new char  OUT_FILE_STR.c_str();
        VARIANCE = 2.0 * (float) dimensionality * D;
        clog << "database : " << OUT_FILE_STR << endl;
        clog << "D  : " << D << endl;
        clog << "k  : " << k << endl;
        clog << "dt : " << dt << endl; 
        clog << "variance : " << VARIANCE << endl; 
    }
    ~Params(){
    }

void register_table(sqlite3pp::database * db) const {
    bool drop = false;
    if (drop){
        std::string sql_drop( "DROP TABLE IF EXISTS register");
        db->execute(sql_drop.c_str());
    }

    char sql[] = "CREATE TABLE IF NOT EXISTS register ("  \
              "id  TEXT KEY  , " \
              "particles INT, " \
              "grid_points INT, " \
              "time_points INT, " \
              "diffusion FLOAT , " \
              "young FLOAT, " \
              "dt FLOAT" \
              ");";
    try {
         clog << "[sqlite3 :] " << sql << " >>> " ;
         clog << db->execute(sql) << endl;
    } catch (exception& ex) {
         cerr << ex.what() << endl;
    }
    
    char sql_ins[] = "INSERT OR REPLACE INTO register "  \
                  "(id, particles, grid_points, time_points, diffusion, young, dt) " \
           "VALUES (:id, :particles, :grid_points, :time_points, :diffusion, :young, :dt); ";

    sqlite3pp::command cmd( *db, sql_ins);
    cmd.bind(":id", RUN_NAME_STR.c_str() );
    cmd.bind(":particles", N );
    cmd.bind(":grid_points", (int) GRID_X_POINTS );
    cmd.bind(":time_points", T );
    cmd.bind(":diffusion", D );
    cmd.bind(":young", k );
    cmd.bind(":dt", dt );
    
    try {
         clog << "[sqlite3 :] " << sql_ins << endl ;
         if ( cmd.execute() != 0) {
             clog << RUN_NAME_STR.c_str() << "\t" ; 
             cerr << db->error_msg() << endl;
         }
    } catch (exception& ex) {
         cerr << ex.what() << endl;
    }
}

  };




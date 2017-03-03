#ifndef DYNEARTHSOL3D_MARKERSET_HPP
#define DYNEARTHSOL3D_MARKERSET_HPP

#include <string>


// forward declaration
class BinaryOutput;
class BinaryInput;

class MarkerSet
{

public:
    MarkerSet( const std::string& name );
    MarkerSet( const Param& param, Variables& var, const std::string& name );
    MarkerSet( const Param& param, Variables& var, BinaryInput& bin, const std::string& name );
    ~MarkerSet()
    { 
        delete _id;
        delete _eta; 
        delete _elem;
        delete _mattype;
    }

    static void random_eta( double* ); // class method

    void append_random_marker_in_elem( int el, int mt );
    void append_marker( const double *eta, int el, int mt );
    void remove_marker(int i);
    void resize(const int);
    void write_chkpt_file(BinaryOutput &bin) const;
    void read_chkpt_file(Variables &var, BinaryInput &bin);
    void write_save_file(const Variables &var, BinaryOutput &bin) const;

    inline int get_nmarkers() const { return _nmarkers; }
    inline void set_nmarkers(int n) { _nmarkers = n; }

    inline int get_id(int m) const { return (*_id)[m]; }
    inline void set_id(const int m, const int i) { (*_id)[m] = i; }

    inline int get_elem(int m) const { return (*_elem)[m]; }
    inline void set_elem(const int m, const int e) { (*_elem)[m] = e; }

    inline int get_mattype(int m) const { return (*_mattype)[m]; }
    inline void set_mattype(const int m, const int mt) { (*_mattype)[m] = mt; }

    inline const double *get_eta(int m) const { return (*_eta)[m]; }
    inline void set_eta( const int i, const double r[NDIMS] );

private:
    const std::string _name;

    // Didn't create a data type for an individual marker to follow the "structure of arrays" concept.
    // Number of markers (may change during remeshing)

    int _nmarkers;
    int _reserved_space;
    int _last_id;

    // Barycentric (local) coordinate within the reference element
    shapefn *_eta;
    // Containing element
    int_vec *_elem;
    // Material type
    int_vec *_mattype;
    // Unique id
    int_vec *_id;

    void random_markers( const Param&, Variables& );
    void regularly_spaced_markers( const Param&, Variables& );
    void allocate_markerdata( const int );

    int initial_mattype( const Param&, const Variables&,
                         int elem, const double eta[NODES_PER_ELEM],
                         const double *x=NULL );
    int layered_initial_mattype( const Param& param, const Variables &var,
                                 int elem, const double eta[NODES_PER_ELEM],
                                 const double *x);
    int custom_initial_mattype( const Param& param, const Variables &var,
                                int elem, const double eta[NODES_PER_ELEM],
                                const double *x );

};

void remap_markers(const Param&, Variables &, 
                   const array_t &, const conn_t &);
void advect_hydrous_markers(const Param &, const Variables &, double,
                            MarkerSet &, Array2D<int,1> &);

#endif

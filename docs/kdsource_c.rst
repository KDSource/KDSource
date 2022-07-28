kdsource C API
===================

C library for particle generation with KDSource objects.

The KDSource C API consists in the ``kdsource`` shared library. It provides the required tools for loading previously created KDSource XML files and generating new particles. With it, advanced users can develop their own plug-ins or converters for using the generated particles. The ``kdtool resample`` application, the ``KDSource.comp`` McStas component and the ``KDSource.c`` template for TRIPOLI-4 external source are examples of usage of this library.

The ``kdsource`` user interface is defined in the following four header files:
    * ``kdsource.h``: The following two structures, and their user interfaces, are defined:
        - ``KDSource``: Main KDSource object. Models a KDE particle source.
        - ``MultiSource``: Models a set of overlapped KDSource objects.
    * ``plist.h``: The ``PList`` structure, and its user interface, is defined. ``PList`` is a wrapper for MCPL particle lists.
    * ``geom.h``: The ``Geometry`` and ``Metric`` structures, and their user interfaces, are defined. They define the variable treatment.
    * ``utils.h``: General utilities.

In the ``kdsource.h`` header all other headers are included, so it is only necessary to include the ``kdsource.h`` file to use the ``kdsource`` library. The following ``bash`` commands must be used to compile a C program using the ``kdsource`` library::

    KDSOURCE=/path/to/kdsourceinstall
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$KDSOURCE/lib
    gcc example.c -lkdsource -lmcpl -lm -I$KDSOURCE/include -L$KDSOURCE/lib

Where ``/path/to/kdsourceinstall`` is the path to the directory where the ``KDSource`` package was installed, and ``example.c`` is the program to compile.

``KDSource`` and ``MultiSource`` structures
-------------------------------------------

The ``KDSource`` structure models a KDE Monte Carlo particle source. It is usually created by loading a XML file generated with the Python API. The main functionality of this structure is particle sampling.

The ``MultiSource`` structure models a set of overlapped ``KDSource``'s. In each particle sampling a source is randomly chosen, respecting the relative intensities of the sources, and the particle is sampled from that source. It is also possible to implement source biasing, by fixing the source sampling probability independently from its intensity, and properly adjusting the sampled particle weights.

The definitions, data structures and functions declared in the ``kdsource.h`` header file are described as follows:

::

    #define MAX_RESAMPLES 1000
    #define NAME_MAX_LEN 256

``MAX_RESAMPLES`` defines the maximum number of sampling attempts, meaning that if a valid particle is not obtained afer ``MAX_RESAMPLES`` attempts the execution is terminated. ``NAME_MAX_LEN`` defines the maximum length for a file name to be read within the ``KDSource`` package.

::

    typedef double (*WeightFun)(const mcpl_particle_t* part);

Definition of a weighting function. It can be used within some functionalities to apply weighting or biasing based on a particle parameters.

::

    typedef struct KDSource{
        double J;       // Total current [1/s]
        PList* plist;   // Particle list
        Geometry* geom; // Geometry defining variable treatment
    } KDSource;

``KDSource`` structure. Models a KDE source, which is composed of a particle list (``PList``) and a chosen variable treatment (``Geometry``). It also includes a total current value, in [1/s] units, which is used to calculate relative intensities in a ``MultiSource`` structure.

::

    KDSource* KDS_create(double J, PList* plist, Geometry* geom);

Create ``KDSource`` structure based on ``PList`` and ``Geometry`` instances, previously created. A total current value ([1/s] units) must be defined. Anyway, since it is only used to define relative intensities in ``MultiSource`` objects, a 1 value can be used if multiple sources will not be used.

::

    KDSource* KDS_open(const char* xmlfilename);

Load ``KDSource`` structure from a XML KDSource file named ``xmlfilename``. This file usually is created with the Python API.

::

    int KDS_sample2(KDSource* kds, mcpl_particle_t* part, int perturb, double w_crit, WeightFun bias, int loop);

Main function for particle sampling with a ``KDSource`` structure. The sampled particle is stored in ``part``. It includes the following arguments to configure sampling options:
* ``perturb``: If it is 0, the particles are sampled directly from the MCPL file without modification. Else, a perturbation is applied following the kernel distribution and the corresponding bandwidth, in concordance with the KDE sampling technique.
* ``w_crit``: If it is equal or less than 0, the sampled statistical weight is set as ``w=w_0``, being ``w_0`` the weight of the original particle from the MCPL file. If it is greater than 0, ``w`` is normalized to 1, using the following technique:

    * If ``w_0<w_crit``: The value ``w_0/w_crit`` is used as probability of using the taken particle, instead of skipping it and taking the next one in the list.
    * If ``w_0>w_crit``: The value ``w_crit/w_0`` is used as probability of stepping forward in the list after sampling.

    This way the sampled particle has always ``w=1`` while respecting the particles weights. The recommended ``w_crit`` is the mean weight of the particles in the list.

* ``bias``: Weighting function for sampling biasing. Will be ignored if ``w_crit<=0``.
* ``loop``: If it is 0, ``exit(EXIT_SUCCESS)`` is called when the end of the list is reached. Else, the list is rewound when its end is reached.

::

    int KDS_sample(KDSource* kds, mcpl_particle_t* part);

Easy function for sampling particles with a ``KDSource`` structure. It redirects to ``KDS_sample2``, with arguments ``perturb=1``, ``w_crit=1``, ``bias=NULL`` and ``loop=1``.

::

    double KDS_w_mean(KDSource* kds, int N, WeightFun bias);

Compute mean particle weight in the list used by the source pointed by ``kds``. ``N`` particles are used in the calculation. If a non-``NULL`` ``bias`` function is set, a bias function can be included.

::

    void KDS_destroy(KDSource* kds);

Destroy ``KDSource`` structure, freeing all associated memory.

::

    typedef struct MultiSource{
        int len;      // Number of sources
        KDSource** s; // Array of sources
        double J;     // Total current [1/s]
        double* ws;   // Frequency weights of sources
        double* cdf;  // cdf of sources weights
    } MultiSource;

``MultiSource`` structure. Models a set of overlapped KDE sources. The values in the ``ws`` array define the sampling frequencies for each source, while the intensities are obtained from the ``J`` parameter of each ``KDSource``.

::

    MultiSource* MS_create(int len, KDSource** s, const double* ws);

Create ``MultiSource`` structure based on the amount of sources, the array of ``KDSource`` structures, and the desired sampling frequencies. During creation the total current ``J`` is computed, as well as the cumulative density function ``cdf``.

::

    MultiSource* MS_open(int len, const char** xmlfilenames, const double* ws);

Load set of ``KDSource`` structures from the KDSource XML files defined in ``xmlfilenames``, and build ``MultiSource`` structure.

::

    int MS_sample2(MultiSource* ms, mcpl_particle_t* part, int perturb, double w_crit, WeightFun bias, int loop);

Main function for particle sampling with a ``MultiSource`` structure. A ``KDSource`` structure is randomly chosen using the frequencies defined in ``ms->ws``, and its sampling function is called passing it the same arguments. After the sampling the obtained particle weight is multiplied by the following source-biasing factor:

.. math::

    f_{SB} = \frac{J_i / J_{tot}}{w_i / w_{tot}}

Where the sub-index :math:`i` represents the chosen source, and :math:`tot` means sum over all sources. This way the discrepancy between the relative intensity and relative sampling frequency is corrected, following the source-biasing technique.

::

    int MS_sample(MultiSource* ms, mcpl_particle_t* part);

Easy function for sampling particles with a ``MultiSource`` structure. It redirects to ``MS_sample2``, with arguments ``perturb=1``, ``w_crit=1``, ``bias=NULL`` and ``loop=1``.

::

    double MS_w_mean(MultiSource* ms, int N, WeightFun bias);

Compute mean weight of the particles in the lists of all sources. Computes the mean weight of each source through the ``KDS_w_mean`` function, with the same ``N`` and ``bias`` parameters, and computes the global mean weight as the weighted average using the ``ms->ws`` values as weights.

::

    void MS_destroy(MultiSource* ms);

Destroy ``MultiSource`` structure, freeing all associated memory.


``PList`` structure
-------------------

The ``PList`` structure models a particle list, and acts as a wrapper for a MCPL file. It allows access to the particles, and includes the possibility of applying a translation and a rotation after loading each particle.

The structures and functions declared in the ``plist.h`` file are described as follows.

::

    typedef struct PList{
        char pt;                     // Particle type ("n", "p", "e", ...)
        int pdgcode;                 // PDG code for particle type

        char* filename;              // Name of MCPL file
        mcpl_file_t file;            // MCPL file

        double* trasl;               // PList translation
        double* rot;                 // PList rotation
        int x2z;                     // If true, apply permutation x,y,z -> y,z,x

        const mcpl_particle_t* part; // Pointer to selected particle
    } PList;

Definition of the ``PList`` structure. The particle type is fixed by the ``pt`` parameter. It includes the MCPL file structure allowing access to it, and (optionally), parameters that define a spatial transformation to be applied to each particle after reading it, usually to adapt the reference system between two simulations. The ``part`` parameter always points to the last particle read.

::

    PList* PList_create(char pt, const char* filename, const double* trasl, const double* rot, int switch_x2z);

Create ``PList`` structure. The particle type must be set (``"n"`` for reutron, ``"p"`` for photon, ``"e"`` for electron), as well as the MCPL file name. Optionally, a translation vector ``trasl`` (3D array) and a rotation vector ``rot`` (3D array, axis-angle format) can be defined, or a ``NULL`` value can be used otherwise. Rotation is applied after translation. If ``switch_x2z`` is non-zero, after applying translation and rotation (if any), the following transformation is applied:

.. math::

    (x,y,z) -> (y,z,x)

::

    int PList_get(const PList* plist, mcpl_particle_t* part);

Get particle, apply spatial transformations (if any), and save it in ``part``. The particle is copied from ``plist->part``, and such variable is not modified afterwards.

::

    int PList_next(PList* plist, int loop);

Step forward in the list, updating the ``plist->part`` variable, until the next valid particle. A particle is considered valid if it has a statistical weight greater than zero and its PDG code (particle type) matches the ``PList`` particle type.

::

    void PList_destroy(PList* plist);

Destroy ``PList`` structure, freeing all associated memory.


``Geometry`` and ``Metric`` structures
--------------------------------------

The main task of the ``Geometry`` structure is to perturb particles following the kernel distribution, as well as the user specified change of variables. This is performed by redirecting the task to the ``Metric`` corresponding to each variables set (e.g. energy, position, direction), which use a perturbation function specific to each variable treatment. ``Geometry`` also manages bandwidths and variables normalization (scaling).

The definitions, structures and functions declared in the ``geom.h`` file are described as follows.

::

    typedef struct Metric Metric;

    typedef int (*PerturbFun)(const Metric* metric, mcpl_particle_t* part,
        double bw);

    struct Metric{
        int dim;            // Dimension
        float* scaling;     // Variables scaling
        PerturbFun perturb; // Perturbation function
        int nps;            // Number of metric parameters
        double* params;     // Metric parameters
    };

Definition of the ``Metric`` structure, along with the perturbation function definition. The main task of ``Metric`` is perturbating a set of variables, which is archived by calling the function pointed by ``perturb``, which in turn uses, apart from the bandwidth received as argument, the scaling factors and the metric parameters (if any). The meaning of these last parameters varies with each metric. For example, they can represent maximum, minimum or reference values for variables, or source shape parameters.

::

    Metric* Metric_create(int dim, const double* scaling, PerturbFun perturb, int nps, const double* params);

Create ``Metric`` structure. The metric dimensionality (``dim``), the scaling factors (``scaling``) and the amount and values of the metric parameters (``nps`` and ``params``) must be provided. Furthermore, a perturb function (``perturb``) is required, which is usually chosen from the available perturb functions (see ``_metric_names`` and ``_metric_perturb`` static variables for available functions).

::

    void Metric_destroy(Metric* metric);

Destroy ``Metric`` structure, freeing all associated memory.

::

    typedef struct Geometry{
        int ord;          // Number of submetrics
        Metric** ms;      // Submetrics
        char* bwfilename; // Bandwidth file name
        FILE* bwfile;     // Bandwidth file
        double bw;        // Normalized bandwidth

        double* trasl;    // Geometry translation
        double* rot;      // Geometry rotation
    } Geometry;

Definition of the ``Geometry`` structure. It contains an ``ord`` number of ``Metric``'s stored in ``ms``. If an adaptive bandwidth was chosen, ``Geometry`` manages the reading of the bandwidths file. In any case the ``bw`` parameter stores the present bandwidth. ``KDSource`` ensures that the value of ``bw`` always corresponds to the ``part`` parameter of the ``PList`` structure. The ``trasl`` and ``rot`` parameters represent the spatial location and orientation of the source.

::

    Geometry* Geom_create(int ord, Metric** metrics, double bw, const char* bwfilename,
        const double* trasl, const double* rot);

Create ``Geometry`` structure. The order ``ord`` (number of metrics), the previously created ``Metric`` structures must be provided, and the location and rotation of the source can be set (``trasl`` and ``rot``, 3D arrays) must be provided. For non-adaptive models, the bandwidth value must be passed as the ``bw`` argument and ``bwfilename=NULL`` must be set. For adaptive bandwidths (list of bandwidths stored in a single precision binary file) the bandwidths file name must be provided as the ``bwfilename`` argument, and ``bw`` is ignored.

::

    int Geom_perturb(const Geometry* geom, mcpl_particle_t* part);

Perturb particle, with a perturbation that follows the kernel distribution, with de bandwidth given by ``geom->bw``.

::

    int Geom_next(Geometry* geom, int loop);

Step forward in the bandwidth list, for the adaptive bandwidth case. If the bandwidth is non-adaptive the function performs no action.

::

    void Geom_destroy(Geometry* geom);

Destroy ``Geometry`` structure, freeing all associated memory.

::

    #define E_MIN 1e-11
    #define E_MAX 20

Maximum and minimum energy values. If after a perturbation an energy value is outside this range, the perturbation is repeated.

::

    int E_perturb(const Metric* metric, mcpl_particle_t* part, double bw);

Perturb energy, with simple energy metric. In this case ``bw`` times the corresponding ``scaling`` element has energy units (MeV).

::

    int Let_perturb(const Metric* metric, mcpl_particle_t* part, double bw);

Perturb energy, with lethargy metric. In this case ``bw`` times the corresponding ``scaling`` element has lethargy units (dimensionless).

::

    int Vol_perturb(const Metric* metric, mcpl_particle_t* part, double bw);

Perturb position in its 3 dimensions, with simple volumetric position metric. In this case ``bw`` times each corresponding ``scaling`` element has position units (cm).

::

    int SurfXY_perturb(const Metric* metric, mcpl_particle_t* part, double bw);

Perturb position in its dimensions ``x`` and ``y``, with simple surface position metric. In this case ``bw`` times each corresponding ``scaling`` element has position units (cm).

::

    int Guide_perturb(const Metric* metric, mcpl_particle_t* part, double bw);

Perturb position and direction, with neutron guide metric. In this case ``bw`` times the first two ``scaling`` elements has position units (cm), while for the las two elements it has angle units (degrees). The metric for direction is polar, relative to each mirror normal direction. Internally it transform the position and direction variables to guide variables :math:`z,t,\theta,\phi`, perturbs these variables, and anti-transforms back.

::

    int Isotrop_perturb(const Metric* metric, mcpl_particle_t* part, double bw);

Perturb direction, with isotropic metric. In this case ``bw`` times the corresponding ``scaling`` element has angle units (degree). The perturbation follows the so-called Von Mises-Fischer distribution, directional equivalent of the gaussian distribution.

::

    int Polar_perturb(const Metric* metric, mcpl_particle_t* part, double bw);

Perturb direction, with polar metric relative to :math:`z`. In this case ``bw`` times each corresponding ``scaling`` element has angle units (degree). Internally it transforms the direction to :math:`\theta,\phi`, perturbs these variables and anti-transforms back.

::

    int PolarMu_perturb(const Metric* metric, mcpl_particle_t* part, double bw);

Perturb direction, with polar metric relative to :math:`z`, using :math:`\mu=cos(\theta)`. In this case ``bw`` times the first ``scaling`` element has :math:`\mu` units (dimensionless), while times the second element it has angle units (degree). Internally it transforms the direction to :math:`\mu,\phi`, perturbs these variables and anti-transforms back.

::

    static const int _n_metrics = 8;
    static const char *_metric_names[] = {"Energy", "Lethargy", "Vol", "SurfXY", "Guide", "Isotrop", "Polar", "PolarMu"};
    static const PerturbFun _metric_perturbs[] = {E_perturb, Let_perturb, Vol_perturb, SurfXY_perturb, Guide_perturb, Isotrop_perturb, Polar_perturb, PolarMu_perturb};

Static variables storing the amount, names and perturbation functions of the implemented metrics.


General utilities
-----------------

Apart from the previously described structures, a set of general utilities for specific problems is provided. It includes mathematical functions, conversion between particle formats and dosimetric factors.

The functions declared in ``utils.h`` are described as follows.

::

    double rand_norm();

Get random value with normal distribution, with zero mean and unit dispersion. Internally it uses the Box-Muller method.

::

    double *traslv(double *vect, const double *trasl, int inverse);
    double *rotv(double *vect, const double *rotvec, int inverse);

Translate and rotate 3D vector. It performs in-place transformation over ``vect`` and returns it. If ``inverse`` is non-zero the inverse transformation is applied. ``trasl`` is the translation vector, while ``rot`` is the rotation vector in axis-angle format.

::

    int pt2pdg(char pt);
    char pdg2pt(int pdgcode);

Convert particle type from ``char`` format (``"n"``: neutron, ``p``: photon, ``e``: electron) to PDG code, and vice versa. An invalid ``pt`` value is converted to the 0 PDG code, and an invalid ``pdgcode`` value is converted to the ``"0"`` character.

::

    double interp(double x, const double *xs, const double *ys, int N);

Interpolation function. The arrays ``xs`` and ``ys``, of length ``N``, store the values to be interpolated. The value interpolated at position ``x`` is returned.

::

    double H10_n_ARN(double E);
    double H10_p_ARN(double E);
    double H10_n_ICRP(double E);
    double H10_p_ICRP(double E);

Dosimetric factors as a function of energy, with units :math:`[pSv cm^2]`. ``n`` indicates neutron and ``p`` indicates photon. ``ARN`` and ``ICRP`` indicates the reference for the interpolation table to be used. The interpolation is performed in logarithmic scale.

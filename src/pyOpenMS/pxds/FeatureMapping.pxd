from Types cimport *
from BaseFeature cimport *
from MSExperiment cimport *
from KDTreeFeatureMaps cimport *

cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>" namespace "OpenMS":

    cdef cppclass FeatureMapping "OpenMS::FeatureMapping":
       FeatureMapping() nogil except +
       FeatureMapping(FeatureMapping) nogil except +

    cdef cppclass FeatureMapping_FetureMappingInfo "OpenMS::FeatureMapping::FetureMappingInfo":
       FeatureMapping_FetureMappingInfo() nogil except +
       FeatureMapping_FetureMappingInfo(FeatureMapping_FetureMappingInfo) nogil except +

    cdef cppclass FeatureMapping_FeatureToMs2Indices "OpenMS::FeatureMapping::FeatureToMs2Indices":
        FeatureMapping_FeatureToMs2Indices() nogil except +
        FeatureMapping_FeatureToMs2Indices(FeatureMapping_FeatureToMs2Indices) nogil except +
        
cdef extern from "<OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>" namespace "OpenMS::FeatureMapping":

       # wrap static method:
       FeatureMapping_FeatureToMs2Indices assignMS2IndexToFeature(MSExperiment& spectra,
                                                                  KDTreeFeatureMaps& fp_map_kd,
                                                                  double precursor_mz_tolerance,
                                                                  double precursor_rt_tolerance,
                                                                  bool ppm) nogil except +  # wrap-attach:FeatureMapping

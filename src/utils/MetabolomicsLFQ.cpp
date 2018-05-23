// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka, Fabian Aicheler $
// $Authors: Oliver Alka, Fabian Aicheler $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>

#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>

#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_MetabolomicsLFQ
 **/

class UTILMetabolomicsLFQ :
  public TOPPBase
{
public:
  UTILMetabolomicsLFQ() :
    TOPPBase("MetabolomicsLFQ", "A standard metabolomics LFQ pipeline.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<file list>", StringList(), "input files");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerInputFile_("design", "<file>", "", "design file");
    setValidFormats_("design", ListUtils::create<String>("tsv"));

    registerOutputFile_("out", "<file>", "", "output mzTab file");
    setValidFormats_("out", ListUtils::create<String>("tsv")); // TODO: add file extension for mzTab

    /// TODO: think about export of quality control files (qcML?)

    Param pp_defaults = PeakPickerHiRes().getDefaults();
//    Param ma_defaults = MapAlignmentPoseClustering().getDefaults();
//    Param fl_defaults = FeatureGroupingAlgorithmKD().getDefaults();

    Param combined;
    
    // Default parameter for FeatureFinderMetabo 
    Param p_com;
    p_com.setValue("noise_threshold_int", 10.0, "Intensity threshold below which peaks are regarded as noise.");
    p_com.setValue("chrom_peak_snr", 3.0, "Minimum signal-to-noise a mass trace should have.");
    p_com.setValue("chrom_fwhm", 5.0, "Expected chromatographic peak width (in seconds).");
    combined.insert("Quantification_common:", p_com);
    combined.setSectionDescription("Quantification_common", "Common parameters for Mass Trace Detection, Elution Profile Detection and Feature Finding");

    Param p_mtd = MassTraceDetection().getDefaults();
    p_mtd.remove("noise_threshold_int");
    p_mtd.remove("chrom_peak_snr");
    combined.insert("Quantification_mtd:", p_mtd);
    combined.setSectionDescription("Quantification_mtd", "Mass Trace Detection parameters");

    Param p_epd;
    p_epd.setValue("enabled", "true", "Enable splitting of isobaric mass traces by chromatographic peak detection. Disable for direct injection.");
    p_epd.setValidStrings("enabled", ListUtils::create<String>("true,false"));
    p_epd.insert("", ElutionPeakDetection().getDefaults());
    p_epd.remove("chrom_peak_snr");
    p_epd.remove("chrom_fwhm");

    combined.insert("Quantification_epd:", p_epd);
    combined.setSectionDescription("Quantification_epd", "Elution Profile Detection parameters (to separate isobaric Mass Traces by elution time).");

    Param p_ffm = FeatureFindingMetabo().getDefaults();
    p_ffm.remove("chrom_fwhm");
    p_ffm.remove("report_chromatograms");
    combined.insert("Quantification_ffm:", p_ffm);
    combined.setSectionDescription("Quantification_ffm", "FeatureFinder parameters (assembling mass traces to charged features)");

    combined.insert("Centroiding:", pp_defaults);
    //combined.insert("Alignment", ma_defaults);
    //combined.insert("Linking:", fl_defaults);

    registerFullParam_(combined);
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parameter handling
    //-------------------------------------------------------------

    // TODO: check handling of single MS file and n-MS files of a single run

    // Read tool parameters
    StringList in = getStringList_("in");
    String out = getStringOption_("out");
    String design_file = getStringOption_("design");
    
    ExperimentalDesign design = ExperimentalDesign::load(design_file);
    std::map<unsigned int, std::vector<String> > frac2ms = design.getFractionToMSFilesMapping();
    
    // Parameter for PeakPickerHiRes
    Param pp_param = getParam_().copy("Centroiding:", true);
    writeDebug_("Parameters passed to PeakPickerHiRes algorithm", pp_param, 3);
    PeakPickerHiRes pp;
    pp.setLogType(log_type_);
    pp.setParameters(pp_param);
   
    // Parameter for FeatureFinderMetabo
    Param common_param = getParam_().copy("Quantification_common:", true);
    writeDebug_("Common parameters passed to sub-algorithms (mtd and ffm)", common_param, 3);

    Param mtd_param = getParam_().copy("Quantification_mtd:", true);
    writeDebug_("Parameters passed to MassTraceDetection", mtd_param, 3);

    Param epd_param = getParam_().copy("Quantification_epd:", true);
    writeDebug_("Parameters passed to ElutionPeakDetection", epd_param, 3);

    Param ffm_param = getParam_().copy("Quantification_ffm:", true);
    writeDebug_("Parameters passed to FeatureFindingMetabo", ffm_param, 3);

    // Parameter for MapAlignPoseClustering
    /*
    Param ma_param = getParam_().copy("Alignment:", true);
    writeDebug_("Parameters passed to MapAlignmentPoseClustering algorithm", ma_param, 3);
    MapAlignmentPoseClustering ma;
    ma.setLogType(log_type_);
    ma.setParameters(ff_param);
    
    // Parameter for FeautreGrouingAlgorithmKD
    Param fl_param = getParam_().copy("Linking:", true);
    writeDebug_("Parameters passed to FeatureGroupingAlgorithmKD algorithm", fl_param, 3);
    FeatureGroupingAlgorithmKD fl;
    fl.setLogType(log_type_);
    fl.setParameters(ff_param);
    */
    
    //-------------------------------------------------------------
    // Loading input
    //-------------------------------------------------------------

    ConsensusMap consensus;
    for (auto const ms_files : frac2ms) // for each fraction->ms file(s)
    {
      vector<FeatureMap> feature_maps;

      //TODO: check if we want to parallelize that
      for (String const & mz_file : ms_files.second) // for each MS file
      {
        // TODO: check if s is part of in 
      
        // load raw file
        MzMLFile mzML_file;
        mzML_file.setLogType(log_type_);

        PeakMap ms_raw;
        mzML_file.load(mz_file, ms_raw);

        if (ms_raw.empty())
        {
          LOG_WARN << "The given file does not contain any spectra.";
          return INCOMPATIBLE_INPUT_DATA;
        }

        // check if spectra are sorted
        for (Size i = 0; i < ms_raw.size(); ++i)
        {
          if (!ms_raw[i].isSorted())
          {
            ms_raw[i].sortByPosition();
            writeLog_("Info: Sorte peaks by m/z.");
          }
        }

        //-------------------------------------------------------------
        // Centroiding of MS1
        //-------------------------------------------------------------
        // TODO: only pick if not already picked (add auto mode that skips already picked ones)
        PeakMap ms_centroided;
        pp.pickExperiment(ms_raw, ms_centroided, true);
        ms_raw.clear(true);// free memory of profile PeakMaps

        // writing picked mzML files for data submission
        // annotate output with data processing info
        // TODO: how to store picked files? by specifying a folder? or by output files that match in number to input files
        // TODO: overwrite primaryMSRun with picked mzML name (for submission)
        // mzML_file.store(OUTPUTFILENAME, ms_centroided);
        // TODO: free all MS2 spectra (to release memory!)

        //-------------------------------------------------------------
        // Feature detection
        //-------------------------------------------------------------
      
        ms_centroided.sortSpectra(true);
        vector<MassTrace> m_traces;
     
        // FeatureFinderMetabo: mass trace detection 
        MassTraceDetection mtdet;
        mtd_param.insert("", common_param);
        mtd_param.remove("chrom_fwhm");
        mtdet.setParameters(mtd_param);

        mtdet.run(ms_centroided, m_traces);

        // FeatureFinderMetabo: elution peak detection
        std::vector<MassTrace> m_traces_final;
        if (epd_param.getValue("enabled").toBool())
        {   
          std::vector<MassTrace> splitted_mtraces;
          epd_param.remove("enabled"); // artificially added above
          epd_param.insert("", common_param);
          ElutionPeakDetection epdet;
          epdet.setParameters(epd_param);
          // fill mass traces with smoothed data as well .. bad design..
          epdet.detectPeaks(m_traces, splitted_mtraces);
          if (epdet.getParameters().getValue("width_filtering") == "auto")
          {
            m_traces_final.clear();
            epdet.filterByPeakWidth(splitted_mtraces, m_traces_final);
          }
          else
          {
            m_traces_final = splitted_mtraces;
          }
        }
        else // no elution peak detection
        {
          m_traces_final = m_traces;
          for (Size i = 0; i < m_traces_final.size(); ++i) // estimate FWHM, so .getIntensity() can be called later
          {
            m_traces_final[i].estimateFWHM(false);
          }
          if (ffm_param.getValue("use_smoothed_intensities").toBool())
          {
            LOG_WARN << "Without EPD, smoothing is not supported. Setting 'use_smoothed_intensities' to false!" << std::endl;
            ffm_param.setValue("use_smoothed_intensities", "false");
          }
        }
        
        // FeatureFinderMetabo: run feature finding
        ffm_param.insert("", common_param);
        ffm_param.remove("noise_threshold_int");
        ffm_param.remove("chrom_peak_snr");
        // TODO: how to write chromatograms for multiple files
        // String report_chromatograms = out_chrom.empty() ? "false" : "true";
        String report_chromatograms = false;
        ffm_param.setValue("report_chromatograms", report_chromatograms);       
       
        FeatureMap feat_map;
        StringList ms_runs;
        ms_centroided.getPrimaryMSRunPath(ms_runs);
        feat_map.setPrimaryMSRunPath(ms_runs);
        std::vector< std::vector< OpenMS::MSChromatogram > > feat_chromatograms;

        FeatureFindingMetabo ffmet;
        ffmet.setParameters(ffm_param);
        ffmet.run(m_traces_final, feat_map, feat_chromatograms);
      
        Size trace_count(0);
        for (Size i = 0; i < feat_map.size(); ++i)
        {
          OPENMS_PRECONDITION(feat_map[i].metaValueExists("num_of_masstraces"),
              "MetaValue 'num_of_masstraces' missing from FFMetabo output!");
          trace_count += (Size) feat_map[i].getMetaValue("num_of_masstraces");
        }

        LOG_INFO << "-- FF-Metabo stats --\n"
                 << "Input traces:    " << m_traces_final.size() << "\n"
                 << "Output features: " << feat_map.size() << " (total trace count: " << trace_count << ")" << std::endl;

        // TODO: if chromatograms are needed please see FeatureFinderMetabo.cpp 

        // store ionization mode of spectra (useful for post-processing by AccurateMassSearch tool)
        if (feat_map.size() > 0)
        {
          set<IonSource::Polarity> pols;
          for (Size i = 0; i < ms_centroided.size(); ++i)
          {
            pols.insert(ms_centroided[i].getInstrumentSettings().getPolarity());
          }
          // concat to single string
          StringList sl_pols;
          for (set<IonSource::Polarity>::const_iterator it = pols.begin(); it != pols.end(); ++it)
          {
            sl_pols.push_back(String(IonSource::NamesOfPolarity[*it]));
          }
          feat_map[0].setMetaValue("scan_polarity", ListUtils::concatenate(sl_pols, ";"));
        }


      }

      //-------------------------------------------------------------
      // Align all features of this fraction
      //-------------------------------------------------------------
      // TODO: adjust to MapAlignerPoseClustering
    //  vector<TransformationDescription> transformations;
      //TODO: check if we need to set reference
    //  Size reference_index(0);
    //  aligner.align(feature_maps, transformations, reference_index);

      // find model parameters:
   //   Param model_params = aligner.getDefaults().copy("model:", true);
   //   String model_type = model_params.getValue("type");

  //    if (model_type != "none")
  //    {
  //      model_params = model_params.copy(model_type + ":", true);
  //      for (TransformationDescription t : transformations)
  //      {
  //        t.fitModel(model_type, model_params);
  //      }
  //    }

      // Apply transformations
  //    for (Size i = 0; i < feature_maps.size(); ++i)
  //    {
  //      MapAlignmentTransformer::transformRetentionTimes(feature_maps[i],
  //        transformations[i]);
  //    }                                     

      //-------------------------------------------------------------
      // Link all features of this fraction
      //-------------------------------------------------------------
    }
    //-------------------------------------------------------------
    // Export of MzTab file as final output
    //-------------------------------------------------------------
   
    // TODO: adjust to Metabolomics 
    // Annotate quants to protein(groups) for easier export in mzTab

    // Fill MzTab with meta data and quants annotated in identification data structure

    return EXECUTION_OK;
  }
};


int main(int argc, const char ** argv)
{
  UTILMetabolomicsLFQ tool;
  return tool.main(argc, argv);
}

/// @endcond


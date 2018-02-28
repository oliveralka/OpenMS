#ifndef OPENMS_ANALYSIS_ID_MESSAGEPASSERFACTORY_HPP
#define OPENMS_ANALYSIS_ID_MESSAGEPASSERFACTORY_HPP

#include <cmath>
#include <evergreen/Evergreen/evergreen.hpp>
#include <evergreen/Utility/inference_utilities.hpp>

typedef unsigned long int uiint;

template <typename Label>
class MessagePasserFactory {
private:
    const int minInputsPAF = 3;
    double alpha, beta, gamma, p;
    Label offset;

    inline double notConditionalGivenSum(double summ) {
        return std::pow((1.0 - alpha), summ) * (1.0 - beta);
    }

public:
    TableDependency<Label> createProteinFactor(Label id);

    TableDependency<Label> createPeptideEvidenceFactor(Label id, double prob);

    TableDependency<Label> createSumEvidenceFactor(size_t nrParents, Label nId, Label pepId);

    TableDependency<Label> createSumFactor(size_t nrParents, Label nId);

    AdditiveDependency<Label> createPeptideProbabilisticAdderFactor(const std::set<Label> & parentProteinIDs, Label nId);

    PseudoAdditiveDependency<Label> createBFPeptideProbabilisticAdderFactor(const std::set<Label> & parentProteinIDs, Label nId, const std::vector<TableDependency <Label> > & deps);

    MessagePasserFactory<Label>(double alpha, double beta, double gamma, double p);



    /// Works on a vector of protein indices (potentially not consecutive)
    // TODO we could recollect the protIDs from the union of parents.
    void fillVectorsOfMessagePassers(const std::vector<Label> & protIDs,
                                     const std::vector<std::vector<Label>> & parentsOfPeps,
                                     const std::vector<double> & pepEvidences,
                                     InferenceGraphBuilder<Label> & igb);

    void fillVectorsOfMessagePassersBruteForce(const std::vector<Label> & protIDs,
                                     const std::vector<std::vector<Label>> & parentsOfPeps,
                                     const std::vector<double> & pepEvidences,
                                     InferenceGraphBuilder<Label> & igb);

    //const std::vector<std::set<Label>> getPosteriorVariables(const std::vector<uiint> & protIDs);
    //const std::vector<std::vector<Label>> getPosteriorVariablesVectors(const std::vector<uiint> & protIDs);
    //const std::vector<std::set<Label>> getPosteriorVariables(uiint rangeProtIDs);
};

//IMPLEMENTATIONS:

template <typename L>
MessagePasserFactory<L>::MessagePasserFactory(double alpha_, double beta_, double gamma_, double p_) {
  assert(0 < alpha_ && alpha_ < 1);
  assert(0 < beta_ && beta_ < 1);
  assert(0 < gamma_ && gamma_ < 1);
  //Note: smaller than 1 might be possible but is untested right now.
  assert(p_ >= 1);
  alpha = alpha_;
  beta = beta_;
  gamma = gamma_;
  p = p_;
}

template <typename L>
TableDependency<L> MessagePasserFactory<L>::createProteinFactor(L id) {
  double table[] = {1 - gamma, gamma};
  LabeledPMF<L> lpmf({id}, PMF({0L}, Tensor<double>::from_array(table)));
  return TableDependency<L>(lpmf,p);
}

template <typename L>
TableDependency<L> MessagePasserFactory<L>::createPeptideEvidenceFactor(L id, double prob) {
  double table[] = {1 - prob, prob};
  LabeledPMF<L> lpmf({id}, PMF({0L}, Tensor<double>::from_array(table)));
  return TableDependency<L>(lpmf,p);
}


template <typename L>
TableDependency<L> MessagePasserFactory<L>::createSumEvidenceFactor(size_t nrParents, L nId, L pepId) {
  Tensor<double> table({nrParents + 1 , 2});
  for (unsigned long i=0; i <= nrParents; ++i) {
    double notConditional = notConditionalGivenSum(i);
    u_long indexArr[] = {i,0};
    table[indexArr] = notConditional;
    u_long indexArr2[] = {i,1};
    table[indexArr2] = 1.0 - notConditional;
  }
  //std::cout << table << std::endl;
  LabeledPMF<L> lpmf({nId, pepId}, PMF({0L,0L}, table));
  //std::cout << lpmf << std::endl;
  return TableDependency<L>(lpmf,p);
}

template <typename L>
TableDependency<L> MessagePasserFactory<L>::createSumFactor(size_t nrParents, L nId) {
  Tensor<double> table({nrParents+1});
  for (unsigned long i=0; i <= nrParents; ++i) {
    table[i] = 1.0/(nrParents+1);
  }
  //std::cout << table << std::endl;
  LabeledPMF<L> lpmf({nId}, PMF({0L}, table));
  //std::cout << lpmf << std::endl;
  return TableDependency<L>(lpmf,p);
}

template <typename L>
AdditiveDependency<L> MessagePasserFactory<L>::createPeptideProbabilisticAdderFactor(const std::set<L> & parentProteinIDs, L nId) {
  std::vector<std::vector<L>> parents;
  std::transform(parentProteinIDs.begin(), parentProteinIDs.end(), std::back_inserter(parents), [](const L& l){return std::vector<L>{l};});
  return AdditiveDependency<L>(parents, {nId}, p);
}

template <typename L>
PseudoAdditiveDependency<L> MessagePasserFactory<L>::createBFPeptideProbabilisticAdderFactor(const std::set<L> & parentProteinIDs, L nId, const std::vector<TableDependency<L>> & deps) {
  std::vector<std::vector<L>> parents;
  std::transform(parentProteinIDs.begin(), parentProteinIDs.end(), std::back_inserter(parents), [](const L& l){return std::vector<L>{l};});
  return PseudoAdditiveDependency<L>(parents, {nId}, deps, p);
}

/// Works on a vector of protein indices (potentially not consecutive)
// TODO we could recollect the protIDs from the union of parents.
template <typename L>
void MessagePasserFactory<L>::fillVectorsOfMessagePassers(const std::vector<L> & protIDs,
                                                          const std::vector<std::vector<L>> & parentsOfPeps,
                                                          const std::vector<double> & pepEvidences,
                                                          InferenceGraphBuilder<L> & igb)
{
  //TODO asserts could be loosened
  assert(parentsOfPeps.size() == pepEvidences.size());
  for (std::vector<uiint> parents : parentsOfPeps)
    for (L parent : parents)
      assert(std::find(protIDs.begin(), protIDs.end(), parent) != protIDs.end());

  for (uiint pid : protIDs)
    igb.insert_dependency(createProteinFactor(pid));

  for (uiint j = 0; j < parentsOfPeps.size(); j++)
  {
    igb.insert_dependency(createPeptideEvidenceFactor(j,pepEvidences[j]));
    igb.insert_dependency(createSumEvidenceFactor(parentsOfPeps[j],j,j));
    igb.insert_dependency(createPeptideProbabilisticAdderFactor(parentsOfPeps[j],j));
  }
}

template <typename L>
void MessagePasserFactory<L>::fillVectorsOfMessagePassersBruteForce(const std::vector<L> & protIDs,
                                                                    const std::vector<std::vector<L>> & parentsOfPeps,
                                                                    const std::vector<double> & pepEvidences,
                                                                    InferenceGraphBuilder<L> & igb)
{
  assert(parentsOfPeps.size() == pepEvidences.size());
  for (std::vector<uiint> parents : parentsOfPeps)
    for (uiint parent : parents)
      assert(std::find(protIDs.begin(), protIDs.end(), parent) != protIDs.end());

  for (uiint pid : protIDs)
    igb.insert_dependency(createProteinFactor(pid));

  for (uiint j = 0; j < parentsOfPeps.size(); j++)
  {
    std::vector<TableDependency<std::string> > deps;
    auto pepdep = createSumEvidenceFactor(parentsOfPeps[j],j,j);
    auto sumdep = createSumFactor(parentsOfPeps[j],j);
    igb.insert_dependency(createPeptideEvidenceFactor(j,pepEvidences[j]));
    igb.insert_dependency(pepdep);
    deps.push_back(sumdep);
    for (auto parent : parentsOfPeps[j]) {
      deps.push_back(createProteinFactor(parent));
    }

    //igb.insert_dependency(createEmptyPeptideProbabilisticAdderFactor(parentsOfPeps[j],j));
    igb.insert_dependency(createBFPeptideProbabilisticAdderFactor(parentsOfPeps[j],j,deps));
  }
}

/*
template <typename L>
const std::vector<std::set<L>> MessagePasserFactory<L>::getPosteriorVariables(const std::vector<L> & protIDs){
    std::vector<std::set<L>> varSets{};
    for (L protID : protIDs){
        std::set<L> varSet{"Pr" + std::to_string(protID)};
        varSets.push_back(varSet);
    }
    return varSets;
}

template <typename L>
const std::vector<std::vector<std::string>> MessagePasserFactory<L>::getPosteriorVariablesVectors(const std::vector<uiint> & protIDs){
  std::vector<std::vector<std::string>> varVecs{};
  for (uiint protID : protIDs){
    std::vector<std::string> varVec{"Pr" + std::to_string(protID)};
    varVecs.push_back(varVec);
  }
  return varVecs;
}

template <typename L>
const std::vector<std::set<std::string>> MessagePasserFactory<L>::getPosteriorVariables(uiint rangeProtIDs){
    std::vector<std::set<std::string>> varSets{};
    for (uiint i=0; i < rangeProtIDs; ++i){
        std::set<std::string> varSet{"Pr" + std::to_string(i)};
        varSets.push_back(varSet);
    }
    return varSets;
}*/

#endif //OPENMS_ANALYSIS_ID_MESSAGEPASSERFACTORY_HPP

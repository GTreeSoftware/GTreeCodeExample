
#ifndef LIBSVMCLASSIFIER_H
#define LIBSVMCLASSIFIER_H
#include <svm.h>
#include "ngtypes/basetypes.h"

struct LIBSVMClassifier
{
//public:
    LIBSVMClassifier(void);
    ~LIBSVMClassifier(void);
    
    bool CreateSVMProblem(MatXd& foreSamples, MatXd& backSamples);
    void DestroySVMProblem(svm_problem&);

    bool Train();
    bool Predict(MatXd& featureList, MatXd& res);
    bool Predict(MatXd& featureList, MatXd& res, std::vector<double> &probability);

    struct svm_parameter param_;		// set by parse_command_line
    struct svm_problem prob_;		// set by read_problem
    struct svm_model *model_;
};

#endif


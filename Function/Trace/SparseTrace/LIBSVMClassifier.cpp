#include "LIBSVMClassifier.h"

LIBSVMClassifier::LIBSVMClassifier(void)
{
    model_ = 0;
    param_.probability = 1;
    prob_.x = 0;
    prob_.y = 0;
    prob_.l = 0;
}

LIBSVMClassifier::~LIBSVMClassifier(void)
{
    svm_destroy_param(&param_);
    svm_free_and_destroy_model(&model_);
    DestroySVMProblem(prob_);
}

bool LIBSVMClassifier::CreateSVMProblem( MatXd& foreSamples, MatXd& backSamples )
{
    if (foreSamples.cols() == 0 || backSamples.cols() == 0) {
        fprintf(stderr,"ERROR: not enough samples\n");
        return false;
    }
    DestroySVMProblem(prob_);
    prob_.l = int(foreSamples.cols() + backSamples.cols());
    prob_.y = new double[prob_.l];
    prob_.x = new svm_node*[prob_.l];

    int featureNum = int(foreSamples.rows());
    for (int i = 0; i < prob_.l; ++i) {
        prob_.x[i] = new svm_node[featureNum + 1];
        prob_.x[i][featureNum].index = -1;
    }
    for (int i = 0; i < foreSamples.cols(); ++i) {
        prob_.y[i] = 1;
        for(int j = 0; j < featureNum; ++j){
            prob_.x[i][j].index = j + 1;
            prob_.x[i][j].value = foreSamples(j, i);
        }
    }
    for (int i = int(foreSamples.cols()); i < prob_.l; ++i) {
        prob_.y[i] = -1;
        for(int j = 0; j < featureNum; ++j){
            prob_.x[i][j].index = j + 1;
            prob_.x[i][j].value = backSamples(j, i - foreSamples.cols());
        }
    }
    return true;
}

void LIBSVMClassifier::DestroySVMProblem(svm_problem& p)
{
    if (p.y) {
        delete[] p.y;
        p.y = 0;
    }
    if (p.x){
        for (int i = 0; i < p.l; ++i) {
            if (p.x[i]) delete[] p.x[i];
            p.x[i] = 0;
        }
    }
    
    delete[] p.x;
    p.x = 0;
}

bool LIBSVMClassifier::Train()
{
    const char *error_msg;
    error_msg = svm_check_parameter(&prob_, &param_);
    if(error_msg)
    {
        fprintf(stderr,"ERROR: %s\n",error_msg);
        return false;
    }
    if(model_) svm_free_and_destroy_model(&model_);
    model_ = svm_train(&prob_,&param_);

    return true;
}

bool LIBSVMClassifier::Predict( MatXd& featureList, MatXd& res )
{
    if(!model_) return false;
    res = MatXd(1, featureList.cols());
    svm_node *test = new svm_node[featureList.rows() + 1];
    test[featureList.rows()].index = -1;
    for (int i = 0; i < featureList.cols(); ++i) {
        for(int j = 0; j < featureList.rows(); ++j){
            test[j].index = j + 1;
            test[j].value = featureList(j, i);
        }
        res(0, i) = svm_predict(model_, test);
    }
    delete[] test;
    return true;
}

bool LIBSVMClassifier::Predict(MatXd& featureList, MatXd& res, std::vector<double> &probability)
{
    if (!model_) return false;
    probability.clear();
    res = MatXd(1, featureList.cols());
    svm_node *test = new svm_node[featureList.rows() + 1];
    test[featureList.rows()].index = -1;
    double *pro = new double[2]; memset(pro, 0, sizeof(double) * 2);
    probability.clear();
    for (int i = 0; i < featureList.cols(); ++i) {
        for (int j = 0; j < featureList.rows(); ++j){
            test[j].index = j + 1;
            test[j].value = featureList(j, i);
        }
        res(0, i) = svm_predict_probability(model_, test, pro);
        probability.push_back(pro[0]);
    }
    delete[] test;
    return true;
}

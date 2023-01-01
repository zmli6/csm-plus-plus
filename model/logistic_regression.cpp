/*
 * @Author: zmli6 liziming033@gmail.com
 * @Date: 2022-09-24 15:29:50
 * @LastEditors: zmli6 liziming033@gmail.com
 * @LastEditTime: 2022-11-17 21:31:44
 * @FilePath: /CSM/model/logistic_regression.cpp
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#include "logistic_regression.h"
#include <fstream>
#include <algorithm>
#include <random>
#include <iostream>
#include <cmath>
#include "common.h"

namespace ANN{
    int LogisticRegression::init(std::unique_ptr<Database> data, std::unique_ptr<Database> predict, int feature_length, float learning_rate, int epochs,Optimization optim , int batch_size){
        CHECK(data->samples.size() == data->labels.size());
        this->m_ = data->samples.size();//m_ is sample's size
        if(this->m_ < 2){
            std::cout << "logistic regression train samples num is too little: " << this->m_ << std::endl;
            return -1;
        }
        if(learning_rate <= 0){
            std::cout << "learning rate must be greater 0: " << learning_rate << std::endl;
            return -1;
        }
        if(epochs < 1){
            std::cout << "number of epochs cannot be zero or a negative number: " << epochs << std::endl;
            return -1;
        }
        this->alpha_ = learning_rate;
        this->epochs_ = epochs;
        this->feature_length_ = feature_length;
        this->data_ = std::move(data);//train data
        this->predict_ = std::move(predict);// predict data
        this->o_.resize(this->m_);//o_ is the predict value
        this->optim_ = optim;
        this->batch_size_ = batch_size;
        this->sampleTime(true);
        return 0;
    }

    void LogisticRegression::sampleTime(bool sample, std::vector<bool> predictLable){
        if(sample){
            double allCheck = 0;
            double NoCheck = 0;
            double Opt = 0;
            double Wos = 0;
            uint sampleSize = this->data_->samples.size();
            for(int i = 0; i < sampleSize; i++){
                allCheck += this->data_->posTime[i];
                NoCheck += this->data_->negTime[i];
                if(this->data_->labels[i] == 1){
                    Opt += this->data_->posTime[i];
                    Wos += this->data_->negTime[i];
                }
                else{
                    Opt += this->data_->negTime[i];
                    Wos += this->data_->posTime[i];
                }
            }
            std::cout << "sample Part" << std::endl;
            std::cout << "allCheck Time -> " << allCheck << std::endl;
            std::cout << "NoCheck Time -> " << NoCheck << std::endl;
            std::cout << "opt Time -> " << Opt << std::endl;
            std::cout << "Wos Time -> " << Wos << std::endl;
        }
        else{
            double allCheck = 0;
            double NoCheck = 0;
            double Opt = 0;
            double Wos = 0;
            double predictTime = 0;
            uint tp = 0;
            uint fp = 0;
            uint tn = 0;
            uint fn = 0;
            uint sampleSize = this->predict_->labels.size();
            CHECK(sampleSize == predictLable.size());
            for(int i = 0; i < sampleSize; i++){
                allCheck += this->predict_->posTime[i];
                NoCheck += this->predict_->negTime[i];
                if(this->predict_->labels[i] == 1){
                    Opt += this->predict_->posTime[i];
                    Wos += this->predict_->negTime[i];
                }
                else{
                    Opt += this->predict_->negTime[i];
                    Wos += this->predict_->posTime[i];
                }
                if(predictLable[i]){
                    predictTime += this->predict_->posTime[i];
                }
                else{
                    predictTime += this->predict_->negTime[i];
                }
                if(predictLable[i]){//predict true
                    if(this->predict_->labels[i] == 1)//actual true
                    {
                        tp++;
                    }
                    else{
                        fp++;
                    }
                }
                else//predict false
                {
                    if(this->predict_->labels[i] == 1){
                        fn++;
                    }
                    else{
                        tn++;
                    }
                }
            }
            std::cout << "predict Part" << std::endl;
            std::cout << "allCheck Time -> " << allCheck << std::endl;
            std::cout << "NoCheck Time -> " << NoCheck << std::endl;
            std::cout << "opt Time -> " << Opt << std::endl;
            std::cout << "Wos Time -> " << Wos << std::endl;
            std::cout << "predict Time -> " << predictTime << std::endl;
            std::cout << "tp -> " << tp << std::endl;
            std::cout << "fp -> " << fp << std::endl;
            std::cout << "fn -> " << fn << std::endl;
            std::cout << "tn -> " << tn << std::endl;
        }
    }

    int LogisticRegression::train(const std::string& model){
        this->w_.resize(this->feature_length_, 0.);//weight resize
        generator_real_random_number(this->w_.data(), this->feature_length_, -0.01f, 0.01f, true);//weight init
        generator_real_random_number(&this->b_, 1, -0.01f, 0.01f);//threshold init

        if(this->optim_ == Optimization::BGD){
            for (int iter = 0; iter < epochs_; ++iter) {
                calculate_gradient_descent();
                auto cost_value = calculate_cost_function();
                //fprintf(stdout, "epochs: %d, cost function: %f\n", iter, cost_value);
                if (cost_value < error_) break;
            }
        }
        else{
            this->random_shuffle_.resize(this->data_->samples.size(), 0);
            for(int i =0; i < this->data_->samples.size(); ++i){
                this->random_shuffle_[i] = i;
            }
            float cost_value = 0.;
            for(int iter = 0; iter < this->epochs_; ++iter){
                std::default_random_engine generator;
                std::shuffle(this->random_shuffle_.begin(), this->random_shuffle_.end(), generator);
            
                int loop = (this->m_ + this->batch_size_ - 1) / this->batch_size_;
                for (int i = 0; i < loop; ++i) {
                    int start = i * batch_size_;
                    int end = start + batch_size_ > m_ ? m_ : start + batch_size_;
                    calculate_gradient_descent(start, end);//do what

                    //get predict value
                    for (int i = 0; i < m_; ++i)
                        o_[i] = calculate_activation_function(calculate_z(data_->samples[i]));

                    //get cost
                    cost_value = calculate_cost_function();
                    //fprintf(stdout, "epochs: %d, loop: %d, cost function: %f\n", iter, i, cost_value);
                    if (cost_value < error_) break;
                }
                if (cost_value < error_) break;
            }
        }
        CHECK(this->store_model(model) == 0);
        return 0;
    }

    int LogisticRegression::store_model(const std::string& model) const
    {
        std::ofstream file;
        file.open(model.c_str(), std::ios::binary);
        if (!file.is_open()) {
            fprintf(stderr, "open file fail: %s\n", model.c_str());
            return -1;
        }

        int length = this->w_.size();
        file.write((char*)&length, sizeof(length));
        file.write((char*)this->w_.data(), sizeof(float) * this->w_.size());
        file.write((char*)&this->b_, sizeof(float));

        file.close();
        return 0;
    }

    int LogisticRegression::load_model(const std::string& model)
    {
        std::ifstream file;
        file.open(model.c_str(), std::ios::binary);
        if (!file.is_open()) {
            fprintf(stderr, "open file fail: %s\n", model.c_str());
            return -1;
        }

        int length{ 0 };
        file.read((char*)&length, sizeof(length));
        this->w_.resize(length);
        this->feature_length_ = length;
        file.read((char*)this->w_.data(), sizeof(float) * this->w_.size());
        file.read((char*)&this->b_, sizeof(float));

        file.close();
        return 0;
    }

    float LogisticRegression::predict(const float* data, int feature_length) const
    {
        

        float value{0.};
        for (int t = 0; t < feature_length_; ++t) {
            value += data[t] * this->w_[t];
        }
        value += this->b_;

        return (calculate_activation_function(value));
    }

    void LogisticRegression::predict()
    {
        int predictSampleSize = this->predict_->samples.size();
        int correctRecord = 0;
        std::vector<bool> predictLabel;
        for(int i = 0; i < predictSampleSize; i++){
            float value{0.};
            for (int t = 0; t < feature_length_; ++t) {
                value += this->predict_->samples[i][t] * this->w_[t];
            }
            value += this->b_;
            if((value > 0.5 && this->predict_->labels[i] == 1)||(value < 0.5 && this->predict_->labels[i] == 0)){
                correctRecord++;
            }
            if(value > 0.5){
                predictLabel.push_back(true);
            }
            else{
                predictLabel.push_back(false);
            }
        }
        float acc = (correctRecord + 0.0) / predictSampleSize;
        std::cout << "predict acc is " << acc << std::endl;
        sampleTime(false, predictLabel);
    }

    float LogisticRegression::calculate_z(const std::vector<float>& feature) const
    {
        float z{0.};
        for (int i = 0; i < feature_length_; ++i) {
            z += this->w_[i] * feature[i];
        }
        z += b_;

        return z;
    }

    float LogisticRegression::calculate_z2(const std::vector<float>& feature, const std::vector<float>& vw) const
    {
        float z{0.};
        for (int i = 0; i < feature_length_; ++i) {
            z += (this->w_[i] - this->mu_ * vw[i]) * feature[i];
        }
        z += this->b_;

        return z;
    }

    float LogisticRegression::calculate_cost_function() const
    {   float J{0.};
        for (int i = 0; i < m_; ++i)
            J += 1./2*std::pow(this->data_->labels[i] - this->o_[i], 2);
        return J/m_;
    }

    float LogisticRegression::calculate_activation_function(float value) const
    {
        switch (activation_func_) {
            case ActivationFunction::Sigmoid:
            default: // Sigmoid
                return (1. / (1. + std::exp(-value))); // y = 1/(1+exp(-value))
        }
    }

    /**
     * @description: loss function --》 add time cost
     * @return {*}
     */
    float LogisticRegression::calculate_loss_function() const
    {
        switch (loss_func_) {
            case LossFunction::MSE:
            default: // MSE
                float value = 0.;
                for (int i = 0; i < m_; ++i) {
                    //method 1:origin
                    value += 1/2.*std::pow(this->data_->labels[i] - this->o_[i], 2);
                    // //method 2:timeDiff as a weight of origin loss
                    // value += 1/2.*std::pow(this->data_->labels[i] - this->o_[i], 2) * this->data_->timeDiff[i];
                    // //method 3: only the error predict add the loss
                    // if((this->o_[i] < 0.5 && this->data_->labels[i] == 1) || (this->o_[i] > 0.5 && this->data_->labels[i] == 0)){
                    //     //method 3.1 : timeDiff as a weight of origin loss
                    //     value += 1/2.*std::pow(this->data_->labels[i] - this->o_[i], 2) * this->data_->timeDiff[i];
                    //     //method 3.2 : timeDiff as the loss
                    //     value += this->data_->timeDiff[i];
                    //     //method 3.3 : origin loss
                    //     value += 1/2.*std::pow(this->data_->labels[i] - this->o_[i], 2);
                    // }
                }
                return value/this->m_;
        }
    }

    /**
     * @description: loss function --> need to add time cost
     * @return {*}
     */
    float LogisticRegression::calculate_loss_function_derivative() const
    {
        switch (loss_func_) {
            case LossFunction::MSE:
            default: // MSE
                float value = 0.;
                for (int i = 0; i < m_; ++i) {
                    //method 1: origin
                    value += (this->o_[i] - this->data_->labels[i]);
                    //method 2: timeDiff as a weight of origin loss
                    // value += (this->o_[i] - this->data_->labels[i]) * this->data_->timeDiff[i];
                    // //method 3: only the error predict add the loss
                    // if((this->o_[i] < 0.5 && this->data_->labels[i] == 1) || (this->o_[i] > 0.5 && this->data_->labels[i] == 0)){
                    //     //method 3.1 : timeDiff as a weight of origin loss
                    //     value += (this->o_[i] - this->data_->labels[i]) * this->data_->timeDiff[i];
                    //     //method 3.2 : timeDiff as the loss
                    //     value += this->data_->timeDiff[i];
                    //     //method 3.3 : origin loss
                    //     value += (this->o_[i] - this->data_->labels[i]);
                    // }
                }
                return value/this->m_;
        }
    }

    /**
     * @description: loss function --> need to add time cost
     * @param {float} predictive_value
     * @param {float} true_value
     * @return {*}
     */
    float LogisticRegression::calculate_loss_function_derivative(float predictive_value, float true_label, int sampleID) const
    {
        switch (loss_func_) {
            case LossFunction::MSE:
            default: // MSE
                float loss = 0.0;
                //method 1: origin
                loss = (predictive_value - true_label) * this->data_->timeDiff[sampleID];
                //loss = -(true_label * log(predictive_value) + (1.0 - true_label) * log(1.0 - predictive_value)) * this->data_->timeDiff[sampleID];
                //method 2: timeDiff as a weight of origin loss
                // result = (predictive_value - true_value) * this->data_->timeDiff[sampleID];
                // //method 3: only the error predict add the loss
                // if((predictive_value < 0.5 && true_value == 1) || (predictive_value > 0.5 && true_value == 0)){
                //     //method 3.1 : timeDiff as a weight of origin loss
                //     result = (predictive_value - true_value) * this->data_->timeDiff[sampleID];
                //     //method 3.2 : timeDiff as the loss
                //     result = this->data_->timeDiff[sampleID];
                //     //method 3.3 : origin loss
                //     result = (predictive_value - true_value);
                // }
                return loss;
        }
    }

    void LogisticRegression::calculate_gradient_descent(int start, int end)
    {
        switch (optim_) {
            case Optimization::Nadam: {
                int len = end - start;
                std::vector<float> m(this->feature_length_, 0.), v(this->feature_length_, 0.), mhat(this->feature_length_, 0.), vhat(this->feature_length_, 0.);
                std::vector<float> z(len, 0.), dz(len, 0.);
                float beta1t = 1., beta2t = 1.;
                for (int i = start, x = 0; i < end; ++i, ++x) {
                    z[x] = calculate_z(this->data_->samples[this->random_shuffle_[i]]);
                    dz[x] = calculate_loss_function_derivative(calculate_activation_function(z[x]), this->data_->labels[this->random_shuffle_[i]]);

                    beta1t *= beta1_;
                    beta2t *= beta2_;

                    for (int j = 0; j < feature_length_; ++j) {
                        float dw = data_->samples[random_shuffle_[i]][j] * dz[x];
                        m[j] = beta1_ * m[j] + (1. - beta1_) * dw; // formula 19
                        v[j] = beta2_ * v[j] + (1. - beta2_) * (dw * dw); // formula 19

                        mhat[j] = m[j] / (1. - beta1t); // formula 20
                        vhat[j] = v[j] / (1. - beta2t); // formula 20

                        this->w_[j] = this->w_[j] - this->alpha_ * (beta1_ * mhat[j] + (1. - this->beta1_) * dw / (1. - beta1t)) / (std::sqrt(vhat[j]) + this->eps_); // formula 33
                    }

                    this->b_ -= (this->alpha_ * dz[x]);
                }
            }
                break;
            case Optimization::NAG: {
                int len = end - start;
                std::vector<float> v(this->feature_length_, 0.);
                std::vector<float> z(len, 0), dz(len, 0);
                for (int i = start, x = 0; i < end; ++i, ++x) {
                    z[x] = calculate_z2(this->data_->samples[this->random_shuffle_[i]], v);
                    dz[x] = calculate_loss_function_derivative(calculate_activation_function(z[x]), this->data_->labels[this->random_shuffle_[i]]);

                    for (int j = 0; j < feature_length_; ++j) {
                        float dw = this->data_->samples[this->random_shuffle_[i]][j] * dz[x];
                        v[j] = mu_ * v[j] + this->alpha_ * dw; // formula 5
                        this->w_[j] = this->w_[j] - v[j];
                    }

                    this->b_ -= (this->alpha_ * dz[x]);
                }
            }
                break;
            case Optimization::AdaMax: {
                int len = end - start;
                std::vector<float> m(this->feature_length_, 0.), u(this->feature_length_, 1e-8), mhat(this->feature_length_, 0.);
                std::vector<float> z(len, 0.), dz(len, 0.);
                float beta1t = 1.;
                for (int i = start, x = 0; i < end; ++i, ++x) {
                    z[x] = calculate_z(this->data_->samples[this->random_shuffle_[i]]);
                    dz[x] = calculate_loss_function_derivative(calculate_activation_function(z[x]), this->data_->labels[this->random_shuffle_[i]]);

                    beta1t *= this->beta1_;

                    for (int j = 0; j < this->feature_length_; ++j) {
                        float dw = this->data_->samples[this->random_shuffle_[i]][j] * dz[x];
                        m[j] = this->beta1_ * m[j] + (1. - this->beta1_) * dw; // formula 19
                        u[j] = std::max(this->beta2_ * u[j], std::fabs(dw)); // formula 24

                        mhat[j] = m[j] / (1. - beta1t); // formula 20

                        // Note: need to ensure than u[j] cannot be 0.
                        // (1). u[j] is initialized to 1e-8, or
                        // (2). if u[j] is initialized to 0., then u[j] adjusts to (u[j] + 1e-8)
                        this->w_[j] = this->w_[j] - this->alpha_ * mhat[j] / u[j]; // formula 25
                    }

                    this->b_ -= (this->alpha_ * dz[x]);
                }
            }
                break;
            case Optimization::Adam: {
                int len = end - start;
                std::vector<float> m(this->feature_length_, 0.), v(this->feature_length_, 0.), mhat(this->feature_length_, 0.), vhat(this->feature_length_, 0.);
                std::vector<float> z(len, 0.), dz(len, 0.);
                float beta1t = 1., beta2t = 1.;
                for (int i = start, x = 0; i < end; ++i, ++x) {
                    z[x] = calculate_z(this->data_->samples[this->random_shuffle_[i]]);
                    dz[x] = calculate_loss_function_derivative(calculate_activation_function(z[x]), this->data_->labels[this->random_shuffle_[i]]);

                    beta1t *= beta1_;
                    beta2t *= beta2_;

                    for (int j = 0; j < this->feature_length_; ++j) {
                        float dw = this->data_->samples[this->random_shuffle_[i]][j] * dz[x];
                        m[j] = this->beta1_ * m[j] + (1. - this->beta1_) * dw; // formula 19
                        v[j] = this->beta2_ * v[j] + (1. - this->beta2_) * (dw * dw); // formula 19

                        mhat[j] = m[j] / (1. - beta1t); // formula 20
                        vhat[j] = v[j] / (1. - beta2t); // formula 20

                        this->w_[j] = this->w_[j] - this->alpha_ * mhat[j] / (std::sqrt(vhat[j]) + this->eps_); // formula 21
                    }

                    this->b_ -= (this->alpha_ * dz[x]);
                }
            }
                break;
            case Optimization::Adadelta: {
                int len = end - start;
                std::vector<float> g(this->feature_length_, 0.), p(this->feature_length_, 0.);
                std::vector<float> z(len, 0.), dz(len, 0.);
                for (int i = start, x = 0; i < end; ++i, ++x) {
                    z[x] = calculate_z(this->data_->samples[this->random_shuffle_[i]]);
                    dz[x] = calculate_loss_function_derivative(calculate_activation_function(z[x]), this->data_->labels[this->random_shuffle_[i]]);

                    for (int j = 0; j < this->feature_length_; ++j) {
                        float dw = this->data_->samples[this->random_shuffle_[i]][j] * dz[x];
                        g[j] = this->mu_ * g[j] + (1. - this->mu_) * (dw * dw); // formula 10

                        //float alpha = std::sqrt(p[j] + eps_) / std::sqrt(g[j] + eps_);
                        float change = -std::sqrt(p[j] + this->eps_) / std::sqrt(g[j] + this->eps_) * dw; // formula 17
                        this->w_[j] = this->w_[j] + change;

                        p[j] = this->mu_ * p[j] +  (1. - this->mu_) * (change * change); // formula 15
                    }

                    this->b_ -= (this->eps_ * dz[x]);
                }
            }
                break;
            case Optimization::RMSProp: {
                int len = end - start;
                std::vector<float> g(this->feature_length_, 0.);
                std::vector<float> z(len, 0), dz(len, 0);
                for (int i = start, x = 0; i < end; ++i, ++x) {
                    z[x] = calculate_z(this->data_->samples[this->random_shuffle_[i]]);
                    dz[x] = calculate_loss_function_derivative(calculate_activation_function(z[x]), this->data_->labels[this->random_shuffle_[i]]);

                    for (int j = 0; j < this->feature_length_; ++j) {
                        float dw = this->data_->samples[this->random_shuffle_[i]][j] * dz[x];
                        g[j] = this->mu_ * g[j] + (1. - this->mu_) * (dw * dw); // formula 18
                        this->w_[j] = this->w_[j] - this->alpha_ * dw / std::sqrt(g[j] + this->eps_);
                    }

                    this->b_ -= (this->alpha_ * dz[x]);
                }
            }
                break;
            case Optimization::AdaGrad: {
                int len = end - start;
                std::vector<float> g(this->feature_length_, 0.);
                std::vector<float> z(len, 0), dz(len, 0);
                for (int i = start, x = 0; i < end; ++i, ++x) {
                    z[x] = calculate_z(this->data_->samples[this->random_shuffle_[i]]);
                    dz[x] = calculate_loss_function_derivative(calculate_activation_function(z[x]), this->data_->labels[this->random_shuffle_[i]]);

                    for (int j = 0; j < feature_length_; ++j) {
                        float dw = this->data_->samples[this->random_shuffle_[i]][j] * dz[x];
                        g[j] += dw * dw;
                        this->w_[j] = this->w_[j] - this->alpha_ * dw / std::sqrt(g[j] + this->eps_); // formula 8
                    }

                    this->b_ -= (this->alpha_ * dz[x]);
                }
            }
                break;
            case Optimization::SGD_Momentum: {
                int len = end - start;
                std::vector<float> v(this->feature_length_, 0.);
                std::vector<float> z(len, 0), dz(len, 0);
                for (int i = start, x = 0; i < end; ++i, ++x) {
                    z[x] = calculate_z(this->data_->samples[this->random_shuffle_[i]]);
                    dz[x] = calculate_loss_function_derivative(calculate_activation_function(z[x]), this->data_->labels[this->random_shuffle_[i]]);

                    for (int j = 0; j < this->feature_length_; ++j) {
                        float dw = this->data_->samples[random_shuffle_[i]][j] * dz[x];
                        v[j] = this->mu_ * v[j] + this->alpha_ * dw; // formula 4
                        this->w_[j] = this->w_[j] - v[j];
                    }

                    this->b_ -= (this->alpha_ * dz[x]);
                }
            }
                break;
            case Optimization::SGD:
            case Optimization::MBGD: {
                int len = end - start;
                std::vector<float> z(len, 0), dz(len, 0);
                for (int i = start, x = 0; i < end; ++i, ++x) {
                    z[x] = calculate_z(this->data_->samples[this->random_shuffle_[i]]);
                    dz[x] = calculate_loss_function_derivative(calculate_activation_function(z[x]), this->data_->labels[this->random_shuffle_[i]], this->random_shuffle_[i]);

                    for (int j = 0; j < this->feature_length_; ++j) {
                        float dw = this->data_->samples[this->random_shuffle_[i]][j] * dz[x];
                        this->w_[j] = this->w_[j] - this->alpha_ * dw;
                    }

                    this->b_ -= (this->alpha_ * dz[x]);
                }
            }
                break;
            case Optimization::BGD:
            default: // BGD
                std::vector<float> z(this->m_, 0), dz(this->m_, 0);
                float db = 0.;
                std::vector<float> dw(this->feature_length_, 0.);
                for (int i = 0; i < m_; ++i) {
                    z[i] = calculate_z(this->data_->samples[i]);
                    this->o_[i] = calculate_activation_function(z[i]);
                    dz[i] = calculate_loss_function_derivative(this->o_[i], this->data_->labels[i], i);

                    for (int j = 0; j < this->feature_length_; ++j) {
                        dw[j] += this->data_->samples[i][j] * dz[i]; // dw(i)+=x(i)(j)*dz(i)
                    }
                    db += dz[i]; // db+=dz(i)
                }

                for (int j = 0; j < this->feature_length_; ++j) {
                    dw[j] /= this->m_;
                    this->w_[j] -= this->alpha_ * dw[j];
                }

                this->b_ -= this->alpha_*(db/this->m_);
        }
    }

    bool LogisticRegression::getModelLoad(){
        return this->modelLoad;
    }

    void LogisticRegression::setModelLoad(bool status){
        this->modelLoad = status;
    }
}
#include <cstdlib>
#include <ctime>
#include <vector>
#include <string>
#include <memory>

namespace ANN {

enum class ActivationFunction {
	Sigmoid // logistic sigmoid function
};

enum class LossFunction {
	MSE // Mean Square Error
};

enum class Optimization {
	BGD, // Batch Gradient Descent
	SGD, // Stochastic Gradient Descent
	MBGD, // Mini-batch Gradient Descent
	SGD_Momentum, // SGD with Momentum
	AdaGrad, // Adaptive Gradient
	RMSProp, // Root Mean Square Propagation
	Adadelta, // an adaptive learning rate method
	Adam, // Adaptive Moment Estimation
	AdaMax, // a variant of Adam based on the infinity norm
	NAG, // Nesterov Accelerated Gradient
	Nadam // Nesterov-accelerated Adaptive Moment Estimation
};

inline char* enum_to_string(Optimization optim)
{
	switch(optim) {
		case Optimization::BGD: return "BGD";
		case Optimization::SGD: return "SGD";
		case Optimization::MBGD: return "MBGD";
		case Optimization::SGD_Momentum: return "SGD_Momentum";
		case Optimization::AdaGrad: return "AdaGrad";
		case Optimization::RMSProp: return "RMSProp";
		case Optimization::Adadelta: return "Adadelta";
		case Optimization::Adam: return "Adam";
		case Optimization::AdaMax: return "AdaMax";
		case Optimization::NAG: return "NAG";
		case Optimization::Nadam: return "Nadam";
		default: return "Invalid optim";
	}
}

struct Database {
	Database() = default;
	std::vector<std::vector<float>> samples; // training set
	std::vector<int> labels; // ground truth labels
	std::vector<double> timeDiff; //timeDiff use to calculate loss
	std::vector<double> posTime;
	std::vector<double> negTime;
};

class LogisticRegression { // two categories
public:
	//LogisticRegression(Optimization optim = Optimization::BGD, int batch_size = 1) : optim_(optim), batch_size_(batch_size) {}
	int init(std::unique_ptr<Database> data, std::unique_ptr<Database> predict ,int feature_length, float learning_rate = 0.00001, int epochs = 1000, Optimization optim = Optimization::BGD, int batch_size = 1);
	int train(const std::string& model);
	int load_model(const std::string& model);
	float predict(const float* data, int feature_length) const; // y = 1/(1+exp(-(wx+b)))
	void predict(); // y = 1/(1+exp(-(wx+b)))
	void set_error(float error) { error_ = error; }
	void set_mu(float mu) { mu_ = mu; }
	void set_eps(float eps) { eps_ = eps; }
	void set_beta1(float beta1) { beta1_ = beta1; }
	void set_beta2(float beta2) { beta2_ = beta2; }
	bool getModelLoad();
	void setModelLoad(bool status);
	void sampleTime(bool sample, std::vector<bool> predictLable = {});

private:
	int store_model(const std::string& model) const;
	float calculate_z(const std::vector<float>& feature) const;  // z(i)=w^T*x(i)+b
	float calculate_z2(const std::vector<float>& feature, const std::vector<float>& vw) const;
	float calculate_cost_function() const;
	static int generate_random(int i) { return std::rand()%i; }

	float calculate_activation_function(float value) const;
	float calculate_loss_function() const;
	float calculate_loss_function_derivative() const;
	float calculate_loss_function_derivative(float predictive_value, float true_value, int sampleID = 0) const;
	void calculate_gradient_descent(int start = 0, int end = 0);

	std::unique_ptr<Database> data_; // train data(images, labels)
	std::unique_ptr<Database> predict_; // predict data
	std::vector<int> random_shuffle_; // shuffle the training data at every epoch
	std::vector<float> o_; // predict value
	int epochs_ = 100; // epochs
	int m_ = 0; // train samples num
	int feature_length_ = 0; // weights length
	float alpha_ = 0.00001; // learning rate
	std::vector<float> w_; // weights
	float b_ = 0.; // threshold
	float error_ = 0.001;
	int batch_size_ = 1; // batch size
	float mu_ = 0.9;
	float eps_ = 1e-8;
	float beta1_ = 0.9;
	float beta2_ = 0.999;
	bool modelLoad = false;

	ActivationFunction activation_func_ = ActivationFunction::Sigmoid;
	LossFunction loss_func_ = LossFunction::MSE;
	Optimization optim_ = Optimization::BGD;
}; // class LogisticRegression

} // namespace ANN
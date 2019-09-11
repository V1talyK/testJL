using AutoGrad

function fun23(x)
    return 0.
end

# function sigmoid1(x)
#     return 1. ./ (1. .+ exp.(-x))
# end
#
# function neural_network(W, x)
#     a1 = sigmoid1(x*W[1])
#     return (a1*W[2])[:]
# end
#
# function neural_network_x(x)
#     a1 = sigmoid1(x*W[1])
#     return (a1*W[2])[:]
# end
m = Dense(2,10,σ)
m(input_point)
neural_network_x(input_point')

function Ax(x)
    return x[2] * sin(pi * x[1])
end

function psy_trial(x, net_out)
    return Ax(x) + x[1] * (1 - x[1]) * x[2] * (1 - x[2]) * net_out
end

function loss_function(W, x, y)
    loss_sum = 0.
    for xi in x
      for yi in y

        input_point = [xi, yi]
        net_out = neural_network(W, input_point)[1]
        #net_out_jacobian = jacobian(neural_network_x,input_point)
        #net_out_hessian = jacobian(jacobian(neural_network_x))(input_point)

        psy_t = psy_trial(input_point, net_out)
        psy_t_jacobian = jacobian(psy_trial)(input_point, net_out)
        psy_t_hessian = jacobian(jacobian(psy_trial))(input_point, net_out)
        gradient_of_trial_d2x = psy_t_hessian[1][1]
        gradient_of_trial_d2y = psy_t_hessian[2][2]
        func = fun23(input_point) # right part function
        err_sqr = ((gradient_of_trial_d2x + gradient_of_trial_d2y) - func)^2
        loss_sum += err_sqr
      end
    end
    return loss_sum
end

input_point = [1,1]
neural_network_x(input_point)
gf1(x) = Tracker.gradient(neural_network_x,x)[1]
jf(x) = Tracker.jacobian(neural_network_x,input_point[:]')
jf1(input_point)


dfg = grad(neural_network_x,input_point);

neural_network_x(input_point)
W = [rand(2, 10), rand(10, 1)]
lmb = 0.001

neural_network(W, [1, 1])
x_space = 1:10
y_space = 1:10
loss_function(W, x_space, y_space)

for i in 100
    loss_grad =  grad(loss_function(W, x_space, y_space))

    W[1] = W[1] - lmb * loss_grad[1]
    W[2] = W[2] - lmb * loss_grad[2]
end



fun45(x) = x
grad(fun45)(1,2)
AutoGrad.jacobian(fun45)
jacobian(fun45,[1, 2])
hessian()

function pred = predict_opls(xtrain,ytrain,xtest,num_OPLS_fact)
[model,stats] = opls(xtrain,ytrain,num_OPLS_fact);
[t, t_ortho, Y_pred] = apply_opls_model(xtrain, ytrain, model, xtest);
pred = Y_pred;
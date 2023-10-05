
import torch
import torch.optim as optim
import numpy

def py_scaling_2d(Sigma):
    k = Sigma.shape[1]
    pz_dim = 1
    for i in range(k):
        pz_dim *= 2
    pz = torch.zeros(pz_dim)
    
    for i in range(pz_dim):
        for ki in range(k):
            for kj in range(k):
                pz[i] += Sigma[ki, kj] * ((i & (1 << ki)) > 0) * ((i & (1 << kj)) > 0)
        
        pz[i] = torch.exp(pz[i])
    
    return pz

def py_mu(pz, p):
    mu = torch.zeros((p, p))
    
    for i in range(p):
        for j in range(p):
            for dim in range(pz.shape[0]):
                mu[i, j] += ((1 << i) & dim > 0) * ((1 << j) & dim > 0) * pz[dim]
    
    return mu

def mom_obj_py(Sigma, mu_hat):
    k = Sigma.shape[1]
    pz = py_scaling_2d(Sigma)
    pz /= torch.sum(pz)
    mu_sig = py_mu(pz, k)
    
    squared_norm = torch.sum((mu_sig - mu_hat) ** 2)
    
    return squared_norm

def mom_grad_py(mu_hat, num_iterations=100, learning_rate=0.01):
    
    # convert mu_hat from numpy to torch
    mu_hat = torch.from_numpy(mu_hat).float()

    # Initialize Sigma with all 0s
    Sigma = torch.zeros_like(mu_hat, requires_grad=True)
    
    optimizer = optim.Adam([Sigma], lr=learning_rate)
    
    for i in range(num_iterations):
        optimizer.zero_grad()
        
        loss = mom_obj_py(Sigma, mu_hat)
        loss.backward()
        
        optimizer.step()
        
        if (i+1) % 10 == 0:
            print(f"Iteration {i+1}: Loss = {loss.item()}")
    
    return Sigma.detach().numpy()


#####################################################################################################################
# This script contains the data from string_distance.R and covariates.R scripts.  We grab the information           #
# From the scripts and separate them into modeling components based on the bdm.jl script.  Then,                    #
# we utilize Michael's model from the bdm.jl script.  Finally, we calculate the posterior statistics and visualize  #
# the trace plots.                                                                                                  #
#####################################################################################################################

begin
    using Pkg
    using CSV, DataFrames#, RCall
    using JLD2, LinearAlgebra, XLSX
    using Dates, ProgressBars # progress
    using Random, Statistics, Distributions, SparseArrays, StatsBase
    using MCMCChains, MCMCDiagnosticTools # for MCMC chain management and plotting
    using StatsPlots
    using Distances
end

#cd("/Volumes/Lexar/Likun/Bayesian Dyadic") #for mac
cd("C:/Users/User/OneDrive - University of Missouri/Andrew's Work/Likun/Bayesian Dyadic") #for windows


time_data = CSV.File("./time_data_nodt.csv") |> DataFrame

cov_df = CSV.File("./cov_data_nodt.csv") |> DataFrame 

#cov_df = CSV.File("./cov_data_nodt_tot.csv") |> DataFrame # Total distance to freshwater

cov_df = cov_df[:, 2:end]

y = cov_df[:,1]

X = cov_df[:,2:9]

#insertcols!(X, 1, :Ones => ones(nrow(X)))


X = Matrix(X)


dsij = cov_df[:,10]
#dtij = ones(Int64, 85905)


H = CSV.File("./h_nodt.csv") |> DataFrame
H = H[:,2:end]

H = Matrix(H)


K = CSV.File("./k_nodt.csv") |> DataFrame

K = K[:,2:end] 

K = Matrix(K)


# Assuming time_data is a DataFrame with columns: "Node Long" and "Node Lat"
n = 415

phi_list = 1:10:(4828.003/3) # Range of phi values.  Find the maximum distance.  Figure out the right step size

R = ones(Float64, n, n, length(phi_list))  # Initialize 3D array with ones

# Loop through phi values
for (k, p) in enumerate(phi_list)
    for i in 1:(n-1)
        lon1 = time_data[i, "longitude"]
        lat1 = time_data[i, "latitude"]
        for j in (i+1):n
            lon2 = time_data[j, "longitude"]
            lat2 = time_data[j, "latitude"]

            # Haversine distance in kilometers
            dist_km = haversine((lat1, lon1), (lat2, lon2))

            r_val = exp(-dist_km / p)
            R[i, j, k] = r_val
            R[j, i, k] = r_val  # Symmetric
        end
    end
    R[:, :, k] += I * 0.0001 
end


################################################################################
# Initializing variables 
################################################################################

#weights = vec([10.678 4.259 3.139 2.666 2.397 2.319 2.286 2.242 2.235 2.205]) # from Vagheesh

include("./bdm.jl")

out = mcmc_both(y, X, H, K, phi_list, R, 8000, dsij) #, dtij,


###
### Construct Dictionary
###

out = load("MCMCoutput.jld2")
#out = load("MCMCoutput2.jld2")


outDict = Dict("beta" => out["betaSave"], 
                "theta" => out["thetaSave"], 
                "eta" => out["etaSave"], 
                "s2_y" => out["s2_ySave"], 
                "s2_eta" => out["s2_etaSave"], 
                "s2_theta" => out["s2_thetaSave"], 
                "gamma" => out["gammaSave"], 
                "phi" => out["phiSave"])

# check correlation of Betas

beta_mat = outDict["beta"][:, 6000:end] |> Matrix

cor(beta_mat')


###
### Compute Posterior Means
###

etaMeans = mean(outDict["eta"], dims = 2)
thetaMeans = mean(outDict["theta"], dims = 2)
betaMeans = mean(outDict["beta"], dims = 2)
s2_yMean = mean(outDict["s2_y"])
s2_etaMean = mean(outDict["s2_eta"])
s2_thetaMean = mean(outDict["s2_theta"])
gammaMean = mean(outDict["gamma"], dims = 2)
phiMean = mean(outDict["phi"])


# save a copy of the beta means for potential surface

betas = outDict["beta"]
betas = betas[:, 6999:end]
betaMeans_r = mean(betas, dims = 2) 
beta_lower = mapslices(x -> quantile(x, 0.025), betas; dims=2)
beta_upper = mapslices(x -> quantile(x, 0.975), betas; dims=2)
betaMeans_r = DataFrame(betaMeans_r = vec(betaMeans_r))

CSV.write("./betameans_r.csv", betaMeans_r)


###
### Construct Chains
###

println("Plotting trace plots.")

nBurn = Int(1)
k = 8000

etaChs = Chains(permutedims(outDict["eta"], [2, 1]))
betaChs = Chains(permutedims(outDict["beta"], [2, 1]))
thetaChs = Chains(permutedims(outDict["theta"], [2, 1]))
gammaChs = Chains(permutedims(outDict["gamma"], [2, 1]))
s2_yCh = Chains(outDict["s2_y"])
s2_etaCh = Chains(outDict["s2_eta"])
s2_thetaCh = Chains(outDict["s2_theta"])
phiCh = Chains(outDict["phi"])

## specifically for visualization
etaChsPlot = Chains(permutedims(outDict["eta"][1:10, nBurn:k], [2, 1]), ["η"*string(i) for i in 1:10])
thetaChsPlot = Chains(permutedims(outDict["theta"][1:10, nBurn:k], [2, 1]), ["θ"*string(i) for i in 1:10])
betaChsPlot = Chains(permutedims(outDict["beta"][:, nBurn:k], [2, 1]), ["β"*string(i) for i in 1:(size(outDict["beta"])[1])])
gammaChsPlot = Chains(permutedims(outDict["gamma"][:, nBurn:k], [2, 1]), ["γ"*string(i) for i in 1:1])
s2_yChPlot = Chains(outDict["s2_y"][nBurn:k], ["σ²(y)"])
s2_etaChPlot = Chains(outDict["s2_eta"][nBurn:k], ["σ²(η)"])
s2_thetaChPlot = Chains(outDict["s2_theta"][nBurn:k], ["σ²(θ)"])
phiChPlot = Chains(outDict["phi"][nBurn:k], ["ϕ"])

###
### Plot Traces
###

savefig(Plots.plot(etaChsPlot), "etaChs.pdf")
savefig(Plots.plot(thetaChsPlot), "thetaChs.pdf")
savefig(Plots.plot(betaChsPlot), "betaChs.pdf")
savefig(Plots.plot(s2_yChPlot), "s2_yCh.pdf")
savefig(Plots.plot(s2_etaChPlot), "s2_etaCh.pdf")
savefig(Plots.plot(s2_thetaChPlot), "s2_thetaCh.pdf")
savefig(Plots.plot(gammaChsPlot), "gammaChs.pdf")
savefig(Plots.plot(phiChPlot), "phichs.pdf")

###
### Krieging Data
###

eta = outDict["eta"]  # size (n × nMCMC)
phi = outDict["phi"]  # vector of length nMCMC


orig_points = time_data[:, 5:6] |> Matrix #original points on the scaled scale

krieg = CSV.File("./for_krieg.csv") |> DataFrame

x_new = krieg[:, 4:5] |> Matrix

nMCMC = 8000


function cov_exp(X1::Matrix, X2::Matrix, phi::Float64)
    n1, n2 = size(X1, 1), size(X2, 1)
    C = Matrix{Float64}(undef, n1, n2)
    for i in 1:n1
        for j in 1:n2
            dist = haversine((X1[i,1], X1[i,2]), (X2[j,1], X2[j,2]))
            C[i,j] = exp(-dist / phi)
        end
    end
    return C
end

# Suppose we have posterior eta samples:
# eta has size (415, nphi)
# Each eta[:, k] corresponds to phi_list[k]

m = size(x_new, 1)
n = size(orig_points, 1)
nphi = size(R, 3)

eta_pred = zeros(m, nphi)

for k in 1:nphi
    phi_k = phi_list[k]
    eta_k = eta[:, k]  # your MCMC posterior sample for this phi
    R_orig = R[:, :, k]
    R_neworig = cov_exp(x_new, orig_points, phi_k)
    
    # Kriging using Cholesky for stability
    L = cholesky(R_orig)
    eta_pred[:, k] = R_neworig * (L \ (L' \ eta_k))
end

eta_mean = mean(eta_pred, dims=2)
eta_lower = mapslices(x -> quantile(x, 0.025), eta_pred; dims=2)
eta_upper = mapslices(x -> quantile(x, 0.975), eta_pred; dims=2)


etaDF = DataFrame(
    x = x_new[:,1],
    y = x_new[:, 2],
    z = vec(eta_mean)
)

CSV.write("./etaDF.csv", etaDF)



###
### Likelihood Function
###

# for each row, find the indices of the 1s
i_index = Int[]
j_index = Int[]

for row in eachrow(H)
    inds = findall(x -> x == 1, row)
    push!(i_index, inds[1])
    push!(j_index, inds[2])
end

function compute_log_likelihood(y, X, betaChs, etaChs, thetaChs, s2_yCh, i_index, j_index)
    n_iter = size(betaChs, 1)
    loglik = zeros(n_iter)

    # Convert all Chains to matrices/vectors
    beta_mat = Array(betaChs)
    eta_mat = Array(etaChs)
    theta_mat = Array(thetaChs)
    s2_y_vec = Array(s2_yCh)

    for iter in ProgressBar(1:n_iter)
        β = beta_mat[iter, :]              
        η = eta_mat[iter, :]         
        θ = theta_mat[iter, :]            
        σ2_y = s2_y_vec[iter]              

        μ = X * β .+ η[i_index] .+ θ[i_index] .+ θ[j_index]
        resid = y .- μ

        loglik[iter] = sum(-0.5 * log(2π * σ2_y) .- 0.5 * resid.^2 ./ σ2_y)
    end

    return loglik
end



loglik = compute_log_likelihood(y, X, betaChs, etaChs, thetaChs, s2_yCh, i_index, j_index)
mean_loglik = mean(loglik[12000:end])

plot(loglik)

k = 8 + 787 + 787 + 1 # betas + etas + thetas + sigma^2
N = length(y)

bic = k * log(N) - 2*log(mean_loglik)



# BIC for the difference in distance to fresh water is 17961.28
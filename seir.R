#packages

library(deSolve)
library(ggplot2)
library(devtools)
library(roxygen2)

#variaveis SEIR-SIR

#simula a dinamica SEIR para humanos e SEI para vetores

#'@param N população humana total
#'@param V população vetorial total
#'@param betta taxa de infecção
#'@param gamma taxa de recuperação
#'@param sigma taxa de transgressão dos expostos para os infectados
#'@param u taxa de mortalidade não relacionada a doença
#'@param b frequencia media que um mosquito pica uma pessoa por dia
#'@param bettav taxa de infecção mosquito vetor
#'@param uv taxa de mortalidade do mosquito
#'@param sigmav taxa de progressão dos mosquitos expostos para infectados
#'
#'@export


# FUNCTION SEIR-SEI

seir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {

#Humans

dS <- (N - S)*u - betta*b*(S*Iv/V)
dE <- (betta * b) * (S*Iv/V) - (u + sigma)*E
dI <- (sigma * E) - (u + gamma)*I
dR <- (gamma * I) - (u * R)

#Vectors

dSv <- V - (bettav * b * Sv*I/N) - (uv * Sv)
dEv <- bettav * b * Sv*I/N - (uv + sigmav)*Ev
dIv <- sigmav * Ev - uv * Iv

 return(list(c(dS, dE, dI, dR, dSv, dEv, dIv)))
})
}

#parameters and initial conditions

#value of parameters
parameters <- c(
betta = 0.4,
gamma = 0.2,
sigma = 0.3,
u = 0.01,
b = 0.5,
uv = 0.1,
bettav = 0.6,
sigmav = 0.2,
N = 500,
V = 2000
)

#initial conditions

state <- c(

#humans

S = 499,
E = 0,
I = 1,
R = 0,

#vectors

Sv = 1999,
Ev = 0,
Iv = 1
)

#time

time <- seq(0, 100, by=1)

#solve the EDOs

output <- ode(y = state, times = time, func = seir_model, parms = parameters)

out_df <- as.data.frame(output)

#plot

#Plot SEIR - Humans

 ggplot(out_df, aes(x = time)) +
  geom_line(linewidth = 1.2, aes(y = S, color = "Susceptibles(S)")) +
  geom_line(linewidth = 1.2, aes(y = E, color = "Exposed(E)")) +
  geom_line(linewidth = 1.2, aes(y = I, color = "Infected(I)")) +
  geom_line(linewidth = 1.2, aes(y = R, color = "Recovered(R)")) +
  labs(title = "SEIR - Simulation for Humans",
       x = "Time (days)", y = "Human Population",
       color = "State") +
   scale_color_manual(values = c("Susceptibles(S)" = "green",
                                 "Exposed(E)" = "purple",
                                 "Infected(I)" = "#ff3333",
                                 "Recovered(R)" = "#0099cc"
                                 )) +
  theme_minimal()








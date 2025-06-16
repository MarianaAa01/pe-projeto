#=====================================================================
# Projeto de PE 2024/2025
# Mariana Carvalho (109974)
#
# Código R para as questões 4, 5, 6, 7, 8, 9 e 10
#
#=====================================================================

#=============================
# Código R para o exercício 4
#=============================
# Parâmetros da distribuição Weibull
lambda <- 29
k <- 7
# 1. Cálculo exato de E(X)
E_exato <- lambda * gamma(1 + 1 / k)
# 2. Gerar a amostra e o valor empírico de E(X)
set.seed(1160)
n <- 5000
amostra <- rweibull(n, shape = k, scale = lambda)
E_empirico <- mean(amostra)
# 3. Desvio absoluto arredondado
desvio <- abs(E_exato - E_empirico)
desvio_arredondado <- round(desvio, 4)
# Print do resultado
#cat("Diferença absoluta:", desvio_arredondado, "\n")
cat("Resultado do exercício 4 = ", desvio_arredondado, "\n")
#_________________________________________________________________________________________________________________________________________
#=============================
# Código R para o exercício 5
#=============================
# Semente dada
set.seed(1091)      
# Número de jogadas a simular (row)
nr_lan <- 30000          
# Número de dados por jogada (col)
nr_dice <- 3        
# Simular todos os lançamentos: matriz n x nr_dice
# Usando o matrix gera todos os resultados de uma vez, o que seria mais eficiente do que correr uma função várias vezes
dice_rolls <- matrix(sample(1:6, nr_lan * nr_dice, replace = TRUE), nrow = nr_lan, ncol = nr_dice)
# Soma dos dados por jogada
somas <- rowSums(dice_rolls)
# Frequências relativas das somas 9 e 10
freq_9 <- mean(somas == 9)
freq_10 <- mean(somas == 10)
# Diferença com 4 casas decimais
diferenca <- round(abs(freq_10 - freq_9), 4)
# Print dos resultados
#print(diferenca)
cat("Resultado do exercício 5 = ", diferenca, "\n")
#_________________________________________________________________________________________________________________________________________
#=============================
# Código R para o exercício 6
#=============================
#                               PASSO 1: VALOR EXATO DE Pn
x <- 5.3
n <- 11
valor_exato_pn <- function(n, x) {
  floor_x <- floor(x)
  soma_termos <- 0
  for (k in 0:floor_x) {
    termo <- (-1)^k * choose(n, k) * (x - k)^n
    soma_termos <- soma_termos + termo
  }
  return(soma_termos / factorial(n))
}
pn_exact <- valor_exato_pn(n, x)
#                               PASSO 2: SIMULAÇÃO COM A SEED=3985 
#                                                E 
#                         APROXIMAÇÃO PELO TLC(TEOREMA DO LIMITE CENTRAL)
#============
# COM A SEED
#============
set.seed(3985)
m <- 110
sn_values <- numeric(m)
for (i in 1:m) {
  sample_x <- runif(n, 0, 1)
  sn_values[i] <- sum(sample_x)
}
pn_sim <- mean(sn_values <= x)
#==========
# COM TLC
#==========
mu_sn <- n * 0.5
var_sn <- n * (1/12)
sd_sn <- sqrt(var_sn)
z <- (x - mu_sn) / sd_sn
pn_tlc <- pnorm(z)
#                               PASSO 3: DESVIO ABSOLUTO ENTRE Pn E Pn,tlc
desvio_tlc <- abs(pn_exact - pn_tlc)
#                               PASSO 4: DESVIO ABSOLUTO ENTRE Pn e Pn,sim
desvio_sim <- abs(pn_exact - pn_sim)
#                               PASSO 5: QUOCIENTE ENTRE OS DESVIOS
quociente <- desvio_tlc / desvio_sim
#                               MOSTRAR RESULTADOS COM 4 CASAS DECIMAIS
#cat("1. Valor exato p_n =", round(pn_exact, 4), "\n")
#cat("2a. Aproximação TLC p_n_TLC =", round(pn_tlc, 4), "\n")
#cat("2b. Aproximação por simulação p_n_sim =", round(pn_sim, 4), "\n")
#cat("3. Desvio TLC = ", round(desvio_tlc, 4), "\n")
#cat("3. Desvio da simulação = ", round(desvio_sim, 4), "\n")
#cat("4. Quociente = ", round(quociente, 4), "\n")
cat("Resultado do exercício 6 = ", round(quociente, 4), "\n")
#_________________________________________________________________________________________________________________________________________
#=============================
# Código R para o exercício 7
#=============================
# Dados do enunciado
n <- 12
soma_x <- 61.17
soma_logx <- 19.45
# Médias
media_x <- soma_x / n
media_logx <- soma_logx / n
# Função para encontrar a raiz (equação de verossimilhança para alpha)
veross <- function(alpha) {
  log(alpha) - digamma(alpha) - log(media_x) + media_logx
}
resultado <- uniroot(veross, interval = c(46.6, 63.4))
alpha_hat <- resultado$root
# Estimar lambda
lambda_hat <- alpha_hat / media_x
modo <- (alpha_hat - 1) / lambda_hat
# Mostrar resultado com 2 casas decimais
#round(modo, 2)
cat("Resultado do exercício 7 = ", round(modo, 2), "\n")
#_________________________________________________________________________________________________________________________________________
#=============================
# Código R para o exercício 8
#=============================
# Dados do enunciado
set.seed(1308)
m <- 1300
n <- 10
mu <- 0.5
sigma <- 1
gamma <- 0.91
z <- qnorm(1 - (1 - gamma) / 2)
# Simulação
cobertura <- replicate(m, {
  amostra <- rnorm(n, mean = mu, sd = sigma)
  media <- mean(amostra)
  erro <- z * sigma / sqrt(n)
  li <- media - erro
  ls <- media + erro
  mu >= li && mu <= ls
})
# Proporção de intervalos que contêm o valor verdadeiro
proporcao <- mean(cobertura)
# Quociente
quociente <- proporcao / gamma
# Resultado final com 4 casas decimais
#round(quociente, 4)
cat("Resultado do exercício 8 = ", round(quociente, 4), "\n")
#_________________________________________________________________________________________________________________________________________
#=============================
# Código R para o exercício 9
#=============================
# Dados do enunciado
set.seed(5707)

m <- 900       #numero de amostras
n <- 25        #dimensão da população X
alfa <- 0.08
h0 <- 5
h1 <- 6.2
# Aqui calculamos o valor crítico
# Este é o limiar de decisão (rejeitamos h0 de T0>critico)
c <- qchisq(1 - alfa, df = 2 * n)
# Inicializar o contar de erros_tipo_II
# É incrementado quando não rejeitamos h0 (mesmo que he=True)
contador_erros_II <- 0
#Loop que corre n vezes
for (i in 1:m) {
  #gera uma amostra de tamanho n=25 com média h1=6.2
  amostra <- rexp(n, rate = 1 / h1)  # gerar sob H1
  #calculamos a média amostral usada na estatística
  media_amostral_X <- mean(amostra)
  #calculamos a estatística T0
  T0 <- (2 * n * media_amostral_X) / h0
  #aqui verificamos se não ultrapassa o critico
  #incrementamos o numero de erros_tipo_II caso não ultrapasse
  if (T0 <= c) {
    contador_erros_II <- contador_erros_II + 1
  }
}
#dividir o numero de erros_tipo_II pelo numero total de simulações
beta_chapeu <- contador_erros_II / m
# Calcular valor teórico de beta
beta <- pgamma(c, shape = n, scale = 2 * h1/h0)
# Quociente
quociente <- round(beta_chapeu/beta, 4)
#print(quociente)
cat("Resultado do exercício 9 = ", quociente, "\n")
#_________________________________________________________________________________________________________________________________________
#=============================
# Código R para o exercício 10
#=============================
# 1. Dados da amostra
dados <- c(2.3, 2.7, 5.2, 0.7, 2.9, 0.6, 2.6, 2.2, 3.8, 0.5, 4.9, 5.4, 3.7, 0.4, 4, 3.6, 2, 0.8, 2.5, 2.8, 1.7, 3.3,
           1.5, 0.4, 6.4, 1.5, 6, 2.1, 0.4, 4.6, 3.1, 4.4, 4, 2.1, 5, 3.3, 4.7, 3.4, 4.3, 4.5, 2.3, 0.5, 4.9, 3.5,
           1.8, 1.9, 2.6, 4.3, 4.6, 5.2, 1.6, 2.8, 2.4, 2.8, 1.8, 3.6, 0.8, 5.1, 1.4, 3.2, 1, 6.3, 3.6, 3.6, 1.8,
           0.9, 4.6, 2.5, 5.8, 0.6, 3.3, 3.2, 6.6, 2.6, 2.5, 1.5, 4.1, 1.7, 2.1, 1.5, 0.4, 4.8, 0.4, 1.5, 4.2, 3.3,
           1.2, 8.1, 2.4, 2.8, 2.1, 6.3, 4.2, 1.3, 6, 1.3, 3.7, 2.5, 6.6, 2.7, 1.4, 2, 0.7, 4.3, 3.4, 4.3, 4, 4,
           0.8, 2.3, 2.5, 5.4, 4.3, 0.5, 3.9, 2.2, 3.4, 1.3, 2.4, 4.7, 2, 1.3, 4.4, 2.9, 2.1, 2.5, 1.6, 2.3, 4.4,
           1.9, 1.9, 1.7, 2, 4.2, 3.4, 3.9, 4.3, 1.3, 2.9, 2.2, 5.1, 2.3, 1.9, 2.9, 5.2, 3.4, 2.6, 2.4, 3.2, 1.3,
           3.1, 5.1, 1.4, 4.2, 0.9, 1.3, 2.1, 2.6, 6.2, 1.6, 2.7, 1.7, 2.3, 3.3, 2.8, 1.2, 2.6, 1.5, 2, 2.8, 2.5, 2,
           1.2, 2.2, 2.6, 2.5, 6, 1.9, 3, 3.8, 1.9, 3.2, 3.1, 1.8, 2.6, 1.9, 3.5, 3.7, 1.8, 2.2, 2, 1.3, 2, 1.1,
           2.2, 3.1, 2.9, 1.3, 0.2, 3.9)
set.seed(4497)
n<-140
# Calcular limites das 6 classes equiprováveis sob H0: Rayleigh(σ=2.2)
sigma <- 2.2
k <- 6
indices_sub <- sample(1:length(dados), n, replace = FALSE)
subamostra <- dados[indices_sub]
probs <- seq(0.2, 1.0, by = 0.2)
limites <- sigma * sqrt(-2 * log(1 - probs))
for (i in 2:length(limites)) {
}
# Agrupar observações e contar frequências observadas
freq_obs <- numeric(k)
for (x in subamostra) {
  if (x <= limites[1]) {
    freq_obs[1] <- freq_obs[1] + 1
  } else if (x <= limites[2]) {
    freq_obs[2] <- freq_obs[2] + 1
  } else if (x <= limites[3]) {
    freq_obs[3] <- freq_obs[3] + 1
  } else if (x <= limites[4]) {
    freq_obs[4] <- freq_obs[4] + 1
  } else {
    freq_obs[5] <- freq_obs[5] + 1
  }
}
# Frequências esperadas 
freq_esp <- n / k
# Teste qui-quadrado
chi2_stat <- sum((freq_obs - freq_esp)^2 / freq_esp)
gl <- k - 1
p_valor <- 1 - pchisq(chi2_stat, gl)
# Decisões nos níveis de significância
niveis <- c(0.01, 0.05, 0.10)
decisoes <- character(length(niveis))
for (i in 1:length(niveis)) {
  if (p_valor < niveis[i]) {
    decisoes[i] <- "Rejeitar H0"
  } else {
    decisoes[i] <- "Não rejeitar H0"
  }
}
cat("\ nDecisoes :\n")
cat(" Nivel 1%: ", decisoes [1] , "\n")
cat(" Nivel 5%: ", decisoes [2] , "\n")
cat(" Nivel 10%: ", decisoes [3] , "\n")
# Determinar a resposta correta
if (all(decisoes == "Não rejeitar H0")) {
  resposta <- "c"
  cat("\nResposta: c) Não rejeitar H0 aos n.s. de 1%, 5% e 10%.\n")
} else if (all(decisoes == "Rejeitar H0")) {
  resposta <- "d"
  cat("\nResposta do exercício 10: : d) Rejeitar H0 aos n.s. de 1%, 5% e 10%.\n")
} else if (decisoes[1] == "Não rejeitar H0" && decisoes[2] == "Rejeitar H0" && decisoes[3] == "Rejeitar H0") {
  resposta <- "a"
  cat("\nResposta do exercício 10: : a) Rejeitar H0 aos n.s. de 5% e 10% e não rejeitar H0 ao n.s. de 1%.\n")
} else if (decisoes[1] == "Não rejeitar H0" && decisoes[2] == "Não rejeitar H0" && decisoes[3] == "Rejeitar H0") {
  resposta <- "b"
  cat("\nResposta do exercício 10: : b) Rejeitar H0 ao n.s. de 10% e não rejeitar H0 aos n.s. de 1% e 5%.\n")
} else {
  resposta <- "e"
  cat("\nResposta do exercício 10: : e) Teste é inconclusivo.\n")
}
#_________________________________________________________________________________________________________________________________________
function p = gauss(mu, sigma) %高斯分布
    p= mu + sigma * sqrt(-2.0 * log(rand)) * sin(2.0 * pi * rand);
end
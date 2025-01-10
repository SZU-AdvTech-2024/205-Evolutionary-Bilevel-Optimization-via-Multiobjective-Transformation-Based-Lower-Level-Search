function p = cauchy_g(mu, gamma)%柯西分布
    p= mu + gamma * tan(pi*(rand - 0.5));
end
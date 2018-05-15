## fisher.test versus phyper
# 0 pathway genes selected: p = 1 (both)
matrix(c(0,100,200,700), nrow = 2) %>% fisher.test(alternative = 'greater') %>% extract('p.value')
1 - phyper(q=0,   m=100, n=900, k=200)
1 - phyper(q=0-1, m=100, n=900, k=200)

# 1 pathway gene selected: p = 1 (both)
matrix(c(1,99,199,701), nrow = 2) %>% fisher.test(alternative = 'greater') %>% extract('p.value')
1 - phyper(q=1-1, m=100, n=900, k=200)
phyper(q=1-1, m=100, n=900, k=200, lower.tail = FALSE)

# 2 pathway genes selected: p = 1 (both)
matrix(c(2,98,198,702), nrow = 2) %>% fisher.test(alternative = 'greater') %>% extract('p.value')
1 - phyper(q=2-1, m=100, n=900, k=200)
phyper(q=2-1, m=100, n=900, k=200, lower.tail = FALSE)

# 3 pathway genes selected: p = 1 (both)
matrix(c(3,97,197,703), nrow = 2) %>% fisher.test(alternative = 'greater') %>% extract('p.value')
1 - phyper(q=3-1, m=100, n=900, k=200)
phyper(q=3-1, m=100, n=900, k=200, lower.tail = FALSE)

# 4 pathway genes selected: 0.9999998 (both)
matrix(c(4,96,196,704), nrow = 2) %>% fisher.test(alternative = 'greater') %>% extract('p.value')
1 - phyper(q=4-1, m=100, n=900, k=200)
phyper(q=4-1, m=100, n=900, k=200, lower.tail = FALSE)

# 10 pathway genes selected: 0.9999998 (both)
matrix(c(10,90,190,710), nrow = 2) %>% fisher.test(alternative = 'greater') %>% extract('p.value')
1 - phyper(q=10-1, m=100, n=900, k=200)
phyper(q=10-1, m=100, n=900, k=200, lower.tail = FALSE)

# 25 pathway genes selected
matrix(c(25,75,175,725), nrow = 2) %>% fisher.test(alternative = 'greater') %>% extract('p.value')
1 - phyper(q=25-1, m=100, n=900, k=200)
phyper(q=25-1, m=100, n=900, k=200, lower.tail = FALSE)

# 50 pathway genes selected: both have same value
matrix(c(50,50,150,750), nrow = 2) %>% fisher.test(alternative = 'greater') %>% extract('p.value')
1 - phyper(q=50-1, m=100, n=900, k=200)
phyper(q=50-1, m=100, n=900, k=200, lower.tail = FALSE)

# 80 pathway genes selected: phyper = 0, fisher.test = 1e-43
matrix(c(80,20,120,780), nrow = 2) %>% fisher.test(alternative = 'greater') %>% extract('p.value')
1 - phyper(q=80-1, m=100, n=900, k=200)
phyper(q=80-1, m=100, n=900, k=200, lower.tail = FALSE)

# 90 pathway genes selected: phyper = 0, fisher.test = 1e-59
matrix(c(90,10,110,790), nrow = 2) %>% fisher.test(alternative = 'greater') %>% extract('p.value')
1 - phyper(q=90-1, m=100, n=900, k=200)
phyper(q=90-1, m=100, n=900, k=200, lower.tail = FALSE)

# 100 (all) pathway genes selected: phyper = 0, fisher.test = 1e-81
matrix(c(100,0,100,800), nrow = 2) %>% fisher.test(alternative = 'greater') %>% extract('p.value')
phyper(q=100-1, m=100, n=900, k=200, lower.tail = FALSE)
1 - phyper(q=100-1,   m=100, n=900, k=200)

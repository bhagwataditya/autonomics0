generate_report <- function(){
  dir.create('report', showWarnings = FALSE)
  
}

cat() <- function(line){
  sprintf('\n%s', line) %>% cat()
}

print_preamble <- function(){
  sink('report/preamble.tex')
  '\\documentclass[a4paper]{article}'         %>% cat()
  '\n\\usepackage{booktabs}'                  %>% cat()
  '\n\\usepackage{parskip}'                   %>% cat()
  '\n\\usepackage{float}'                     %>% cat()
  '\n\\usepackage[margin=1.25cm]{geometry}'   %>% cat()
  '\n\\usepackage{pdflscape}'                 %>% cat()
  '\n\\usepackage{longtable}'                 %>% cat()
  '\n\\usepackage{hyperref}'                  %>% cat()
  sink()
}

print_title_page <- function(title){
  sprintf('\\title{%s}', title) %>% cat()
  \author{
    \begin{tabular}{rl}
    Neha Goswami        & Sample Preparation and Data Analysis         \\
    Aditya M. Bhagwat   & Statistics and Bioinformatics           \\
    Hisham Ben Hamidane & Mass Spectromery and Data Processing  \\
    Johannes Graumann   & Coordination
    \end{tabular}
  }
}

print_preamble()
print_title_page(title = 'Studying the orai1 interactome through three knockdowns')



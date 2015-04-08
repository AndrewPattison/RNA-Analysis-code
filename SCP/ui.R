#UI
create_bam_list <- function(){
  result_list <- list('WT-t0-rep1.bam','WT-t0-rep2.bam','WT-t5-rep1.bam', 'WT-t5-rep2.bam', 'WT-t10-rep1.bam', 'WT-t10-rep2.bam', 'WT-t20-rep1.bam', 'WT-t20-rep2.bam',  'RNA1-1-t0-rep1.bam', 'RNA1-1-t0-rep2.bam','RNA1-1-t5-rep1.bam', 'RNA1-1-t5-rep2.bam', 'RNA1-1-t10-rep1.bam', 'RNA1-1-t10-rep2.bam','RNA1-1-t20-rep1.bam', 'RNA1-1-t20-rep2.bam')
  return(result_list)
}

create_gff_list <- function(){
  gff_files <- list.files(getwd(), pattern = "\\.gff$")
  
  return(gff_files)
}


shinyUI(fluidPage(
  titlePanel("Poly(A) Plotter"),
  em(helpText("created by Andrew Pattison, Jack Xu and Paul Harrison for the Beilharz Lab", align = "right")),
  helpText("
           This app analyses the pencentage population against poly-A tail length in a 
           specific gene or peak, from the selected sample files."), br(),
  
  sidebarLayout(
    sidebarPanel(
      h2("Control Panel"),
      helpText("Select the desired parameters, press Submit to apply 
               the changes and view the result."),
      
      checkboxGroupInput("bam_select", label = h4("Bam files selection"),
                         choices = create_bam_list(), 
                         selected = create_bam_list()),
      selectInput("gff_select", label = h4("gff selection"), 
                  choices = create_gff_list(), selected = as.character(create_gff_list()[1])),
      textInput("gene_name", label = h4("gene name or peak number"), 
                value = "enter gene/peak name"),
      numericInput("n_replicates", label = h4("combine by groups of"),
                   value = 2),
      checkboxInput("merge", label = "merge groups", value = T),
      checkboxInput("legend", label = "display legend", value = T), br(),
      submitButton("Submit"),
      
      h3("save"),
      #textInput("save_name", label = h4("save name"), 
      #    value = "enter text..."),
      downloadButton("downloadPlot", label = "Download Plot")
      ),
    
    mainPanel(
      h2("Plot"),
      #textOutput("text1"),
      em(textOutput("text3")),
      plotOutput("plot", height = "600"),
      verbatimTextOutput("text2")
      
    )
  )
))
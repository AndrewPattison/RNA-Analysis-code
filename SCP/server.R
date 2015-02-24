library(Rsamtools)
source("helper.R")
#source("plot_helper.R")


shinyServer(
  function(input, output){
    output$text1 <- renderText({
      c('Order of bam file processing:', create_bam_list(input$bam_select), input$gff_select)
    })
    output$plot <- renderPlot({
      SCP(bam_select = input$bam_select, gff = input$gff_select, 
          name_list = input$gene_name, number_of_replicates = input$n_replicates, 
          combine = input$merge, plot_mean = input$mean, plot_legend = input$legend)  
    })
    output$text2 <- renderPrint({
      number_of_reads()
      
    })
    output$text3 <- renderText({
      if(substr(input$gene_name,1,4) != "peak"){
        c("The gene contains the peak/s:", peak_finder(gff_file = input$gff_select, name_list = input$gene_name))
      }
      else{
        NULL
      }
    })
    output$downloadPlot <- downloadHandler(
      filename = function(){
        paste(input$gene_name, '.eps', sep='')
      },
      content = function(file){
        setEPS(width = 10)
        postscript(file)
        SCP(bam_select = input$bam_select, gff = input$gff_select, 
            name_list = input$gene_name, number_of_replicates = input$n_replicates, 
            combine = input$merge, plot_mean = input$mean, plot_legend = input$legend)
        dev.off()
      }
    )
  }    
)
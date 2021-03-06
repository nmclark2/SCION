#Version 2.0 beta
#July 2020

# Define UI ----
ui <- fluidPage(
  titlePanel("Spatiotemporal Clustering and Inference of Omics Networks (SC-ION)"),
  sidebarLayout(
    sidebarPanel(h3("Files"),
      # Copy the line below to make a file upload manager
      textInput("text", label = "Working directory", value = "c:/path/to/main/dir"),
      fileInput("file", label = "Target Matrix"),
      fileInput("file2", label = "Regulator Matrix"),
      fileInput("file3", label = "Target List"),
      fileInput("file4", label = "Regulator List"),
      radioButtons("radio", label = "Clustering", choices = c("Temporal","Non-Temporal","Upload","None")),
      fileInput("file5", label = "Clustering Matrix"),
      numericInput("num", label = "Clustering Threshold", value = 0.5),
      fileInput("file6", label = "Clusters File"),
      radioButtons("radio2", label = "Hub Connection", choices = c("Yes","No")),
      actionButton("action", label = "Run")
    ),
    mainPanel(h3("Input files"),
              strong("NOTE: Currently SCION only takes comma separated values (.csv) files as input."),
              p(),
              strong("Working directory"),
              p("Name of your working directory which MUST contain all of your files."),
              strong("Target Matrix"),
              p("Data matrix to use for your gene targets (e.g. transcript data). 
                Columns are samples, rows are genes."),
              strong("Regulator Matrix"),
              p("Data matrix to use for your gene regulators (e.g. protein data). 
                Columns are samples, rows are genes."),
              strong("Target List"),
              p("List of genes you consider targets in your network."),
              strong("Regulator Matrix"),
              p("List of genes you consider regulators in your network (e.g. transcription factors)"),
              strong("Clustering"),
              p("Select if you are clustering using Temporal or Non-Temporal data.
                Select Upload to use your own file.
                Select None if you do not wish to cluster."),
              strong("Clustering Matrix"),
              p("Data matrix to use to cluster your gene targets. 
                This can be the same matrix as Target Matrix, or a different matrix of data if desired.
                Columns are samples, rows are genes.
                This is only required if you choose to cluster using Temporal or Non-Temporal data."),
              strong("Clustering Threshold"),
              p("Indicate the clustering threshold to use. The clustering threshold is how similar the
                genes must be to be clustered together. If you pick a threshold close to 1, you will
                obtain many clusters with high within-cluster similarity. If you pick a threshold close to 0,
                you will obtain few clusters with low within-cluster similarity. We recommend starting with 0.5."),
              strong("Clusters File"),
              p("If uploading your own file with pre-determined clusters, select it here.
                The first column must be the names of your genes, and there must be a column called clusters that contains each cluster number.
                For an example of how this file should be formatted, see the clusters.csv file included in the test data.
                This is only required if you choose Upload as your clustering method."),
              strong("Hub Connection"),
              p("Choose Yes if you want to connect the clusters using their hubs (recommended).
                Choose No if you do not want to use hub connection."),
              h3("Output files"),
              strong("All output files are saved to your working directory."),br(),br(),
              p("If you cluster, a folder is created called Cluster Plots which will contain plots of all of your clusters.
              These are good for visualizing if you need to alter your clustering threshold."),
              p("If you cluster, a file called clusters.csv which contains your clustering results."),
              p("If you cluster, a folder is created called Cluster networks which will contain individual networks for each
                of your clusters. These networks are combined to form the final network."),
              p("A text file called FINAL_NETWORK.txt which contains the final network. This can be imported
                into Cytoscape for visualization.")
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  options(shiny.maxRequestSize=100*1024^2) 
  observe({
    if (input$action > 0){
      SCION(target_genes_file=input$file3$datapath,reg_genes_file=input$file4$datapath,
            target_data_file=input$file$datapath,reg_data_file=input$file2$datapath,
            is_clustering=input$radio, clusters_file=input$file6$datapath, hubs=input$radio2,
            clustering_data_file=input$file5$datapath,threshold=input$num,working_dir=input$text)
    }
  })
}


# Run the app ----
shinyApp(ui = ui, server = server)
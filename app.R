#Version 3.0 
#April 2022

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
      radioButtons("radio3", label = "Normalize Edge Weights", choices = c("Yes","No")),
      numericInput('weight', label = "Edge Weight Cutoff", min=0, max=10, value = 0.33, step=0.1),
      radioButtons("radio", label = "Clustering", choices = c("Temporal:DTW","Non-Temporal:ICA","Non-Temporal:k-means","Upload","None")),
      fileInput("file5", label = "Clustering Matrix"),
      numericInput("num", label = "Clustering Threshold", min = 0, max = 5, value = 0.5, step=0.1),
      fileInput("file6", label = "Clusters File"),
      radioButtons("radio2", label = "Hub Connection", choices = c("Yes","No")),
      numericInput("num2", label = "Number of cores (for parallelization)", min = 1, value=1, step=1),
      actionButton("action", label = "Run"),
      actionButton ('screenshot', HTML ('Screenshot this window')),
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
                Columns are samples, rows are genes. If your regulators contain additional information such as a PTM site, this should be denoted using a . (e.g. SOX2.S35) "),
              strong("Target List"),
              p("List of genes you consider targets in your network."),
              strong("Regulator List"),
              p("List of genes you consider regulators in your network (e.g. transcription factors). If your regulators contain additional information such as a PTM site, this should be denoted using a . (e.g. SOX2.S35) "),
              strong("Normalize Edge Weights"),
              p("Option to normalize the edge weights to the [0,1] range before applying the cutoff (default: Yes). If you would rather not normalize the weights, select No. Note that not normalizing the weights will result in a smaller network, as edges with negative weights will be automatically trimmed."),
              strong("Edge Weight Cutoff"),
              p("Cutoff to use for trimming the network. The value must be 0 or greater. A cutoff of 0 keeps all inferred edges. The larger the cutoff, the less edges are kept. Higher cutoff = more confidence in the inferred edges. We recommend using some sort of cutoff except in the case of very small networks. Default: 0.33"),
              strong("Clustering"),
              p("Select if you are clustering using Temporal or Non-Temporal data.
                For Non-Temporal data, there are two choices: ICA and k-means.
                Select Upload to use your own file.
                Select None if you do not wish to cluster."),
              strong("Clustering Matrix"),
              p("Data matrix to use to cluster your gene targets. 
                Columns are samples, rows are genes.
                This is only required if you choose to cluster using Temporal or Non-Temporal data."),
              strong("Clustering Threshold (DTW and ICA only)"),
              p("Indicate the clustering threshold to use. The clustering threshold is how similar the
                genes must be to be clustered together. The threshold works differently depending on the clustering method selected. Default: 0.5."),
                p("For DTW clustering: Range 0 to 1. If you pick a threshold close to 1, you will
                obtain many clusters with high within-cluster similarity. If you pick a threshold close to 0, you will obtain few clusters with low within-cluster similarity."),
                p("For ICA clustering: Range 0 to 5. If you pick a threshold close to 5, you will obtain few clusters with low within-cluster similarity. If you pick a threshold close to 0, you will obtain many clusters with high within-cluster similarity."),
              p("For k-means clustering, the threshold parameter is not used. The number of clusters is selected based on the dimensions of the input matrix. Multiple values are tested, and then the optimal number of clusters is picked based on the silhouette index."),
              strong("Clusters File"),
              p("If uploading your own file with pre-determined clusters, select it here.
                The first column must be the names of your genes, and there must be a column called clusters that contains each cluster number.
                For an example of how this file should be formatted, see the clusters.csv file included in the test data.
                This is only required if you choose Upload as your clustering method."),
              strong("Hub Connection"),
              p("Choose Yes if you want to connect the clusters using their hubs (recommended).
                Choose No if you do not want to use hub connection."),
              strong("Number of cores (for parallelization)"),
              p("Enter the number of cores to use if you would like to parallelize.
                The default is 1 core, which will not run in parallel.
                If you are unsure of the number of cores of your machine, leave this at 1 to disable parallelization, or use detectCores() or other functionality to determine the number of usable cores."),
              h3("Output files"),
              strong("All output files are saved to your working directory."),br(),br(),
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
            target_data_file=input$file$datapath,reg_data_file=input$file2$datapath, weightthreshold=input$weight,normalize_edges=input$radio3,
            is_clustering=input$radio, clusters_file=input$file6$datapath, hubs=input$radio2,
            clustering_data_file=input$file5$datapath,threshold=input$num,num.cores=input$num2,working_dir=input$text)
    }
  })
  #Save a screenshot of the window when the user clicks the button
  observeEvent(input$screenshot, {
    screenshot()
  }
  )
}


# Run the app ----
shinyApp(ui = ui, server = server)
################################################################################
# UI
#
# The Run step lives in the persistent dashboard sidebar (every other tab
# depends on its result, so it stays visible regardless of which tab is
# active), matching protigy-v2's setup-sidebar convention. Network Diagnostics
# and Visualize are navbar tabs in the main body. Help is shown by default so
# it's readable while filling out the sidebar; a completed run auto-switches
# to Network Diagnostics.
################################################################################

#' @keywords internal
SIDEBAR_WIDTH <- 350

#' @keywords internal
app_UI <- function(request) {
  shinydashboardPlus::dashboardPage(
    skin = "blue-light",
    header = shinydashboard::dashboardHeader(title = paste0("SCION v", utils::packageVersion("SCION"))),
    sidebar = shinydashboard::dashboardSidebar(width = SIDEBAR_WIDTH, runSidebarUI("run")),
    body = shinydashboard::dashboardBody(
      shiny::tags$head(shiny::tags$style(shiny::HTML(sprintf(
        ".main-header .logo {
           width: %1$dpx; text-align: left; padding-left: 15px; font-weight: bold;
         }
         .main-header .navbar { margin-left: %1$dpx; }
         .main-sidebar, .main-sidebar .control-label, .main-sidebar label,
         .main-sidebar .help-block, .main-sidebar .radio label, .main-sidebar .checkbox label,
         .main-sidebar h1, .main-sidebar h2, .main-sidebar h3, .main-sidebar h4,
         .main-sidebar h5, .main-sidebar h6 {
           color: #333333;
         }
         .main-sidebar .sidebar {
           height: calc(100vh - 50px);
           overflow-y: auto;
         }
         .main-sidebar .form-group {
           margin-bottom: 6px;
         }
         .main-sidebar .form-group label {
           margin-bottom: 2px;
         }
         .main-sidebar .checkbox, .main-sidebar .radio {
           margin-top: 0; margin-bottom: 6px;
         }",
        SIDEBAR_WIDTH
      )))),
      shinyjs::useShinyjs(),
      shiny::navbarPage(
        title = "",
        id = "navbar-tabs",
        selected = "help",
        shiny::tabPanel("Help", helpTabUI("help"), value = "help", icon = shiny::icon("circle-question")),
        shiny::tabPanel("Network Diagnostics", diagnosticsTabUI("diagnostics"), value = "diagnostics",
                        icon = shiny::icon("magnifying-glass-chart")),
        shiny::tabPanel("Visualize", visualizeTabUI("visualize"), value = "visualize",
                        icon = shiny::icon("circle-nodes"))
      )
    ),
    controlbar = NULL,
    title = "SCION"
  )
}

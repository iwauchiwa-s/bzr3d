library(shiny)
library(plotly)
library(shinyWidgets)

# Define UI
ui <- shinyUI(pageWithSidebar(
  
  #  Application title
  headerPanel("3-D Belousov-Zhabotinsky Reaction Simulator"),
  
  # Sidebar with sliders
  sidebarPanel(
    
    tags$h5("Caution! Selecting large number for grids and steps will take a long integration time. Please start with small numbers."),
    
    actionButton("sim", "Simulate"),
    
    tags$h4("Parameter Settings"),
    
    sliderInput("ngrid", label = h5("Grid number:"), min = 3, 
                max = 20, value = 5),
    
    sliderInput("nstep", label = h5("Integration steps:"), min = 2, 
                max = 100, value = 10),
    
    sliderInput("rlxp", label = h5("Diffusion parameter:"), min = 0, 
                max = 1, value = 0.8),
    
    sliderInput("racp", label = h5("Chemial reaction parameter:"), min = 0, 
                max = 1, value = 1),
    
    br(),
    
  ),
  
  # Show plots
  mainPanel(
    plotlyOutput("graph")
  )
  
))

server <- function(input, output, session) {
  
  mySim <- eventReactive(input$sim, {
    
    # define data length
    n <- input$ngrid
    ni <- n
    nj <- n
    nk <- n
    nn <- ni * nj* nk
    
    # relaxation for diffusion (0 to 1) 0:no diffusion
    rlx <- input$rlxp
    
    # constants
    ca <- input$racp
    cb <- input$racp
    cc <- input$racp
    
    # integration cycle number
    nc <- input$nstep
    
    # define 3-D data and initial random values
    a <- array(runif(nn,0,1), dim = c(ni,nj,nk))
    b <- array(runif(nn,0,1), dim = c(ni,nj,nk))
    c <- array(runif(nn,0,1), dim = c(ni,nj,nk))
    
    # working matrix
    aw <- array(0, dim = c(ni,nj,nk))
    bw <- array(0, dim = c(ni,nj,nk))
    cw <- array(0, dim = c(ni,nj,nk))
    
    # plot 3D for matrix a and prepare plot colors
    p <- c(1:ni)
    q <- c(1:nj)
    r <- c(1:nk)
    df <- expand.grid(p=p,q=q,r=r)
    
    d <- a
    dim(d) <- c(ni*nj*nk,1)
    df$s <- d
    df$t <- 1
    adf <- df
    
    
    for (istep in 1:nc) { # integration steps
      
      for (i in 1:ni) {
        for (j in 1:nj) {
          for (k in 1:nk) {
            # diffusion (local 3x3x3 grid box averaging)
            a1 <- 0
            b1 <- 0
            c1 <- 0
            # 3x3x3 grid summation
            for (x in (i-1):(i+1)) {
              for (y in (j-1):(j+1)) {
                for (z in (k-1):(k+1)) {
                  a1 <- a1 + a[ (( x - 1 ) %% ni ) + 1, (( y - 1 ) %% nj ) + 1, (( z - 1 ) %% nk ) + 1]
                  b1 <- b1 + b[ (( x - 1 ) %% ni ) + 1, (( y - 1 ) %% nj ) + 1, (( z - 1 ) %% nk ) + 1]
                  c1 <- c1 + c[ (( x - 1 ) %% ni ) + 1, (( y - 1 ) %% nj ) + 1, (( z - 1 ) %% nk ) + 1]
                }
              }
            }
            # averaging 
            a1 <- a1 / 27
            b1 <- b1 / 27
            c1 <- c1 / 27
            
            # relaxation effect for diffusion process
            a2 <- a[i,j,k] + ( a1 - a[i,j,k] ) * rlx
            b2 <- b[i,j,k] + ( b1 - b[i,j,k] ) * rlx
            c2 <- c[i,j,k] + ( c1 - c[i,j,k] ) * rlx
            
            # reactions
            a3 <- a2 + a2 * ( ca * b2 - cc * c2 )
            b3 <- b2 + b2 * ( cb * c2 - ca * a2 )
            c3 <- c2 + c2 * ( cc * a2 - cb * b2 )
            # constrain within 0-1
            if ( a3 < 0 ) {
              a3 <- 0
            }
            if ( a3 > 1 ) {
              a3 <- 1
            }
            if ( b3 < 0 ) {
              b3 <- 0
            }
            if ( b3 > 1 ) {
              b3 <- 1
            }
            if ( c3 < 0 ) {
              c3 <- 0
            }
            if ( c3 > 1 ) {
              c3 <- 1
            }
            
            # store change values at grid [i,j,k] into working matrix
            aw[i,j,k] <- a3
            bw[i,j,k] <- b3
            cw[i,j,k] <- c3
          }
        }
      }
      # change matrix data after each step calculation
      a <- aw
      b <- bw
      c <- cw
      
      d <- a
      dim(d) <- c(ni*nj*nk,1)
      df$s <- d
      df$t <- istep + 1
      adf <- rbind(adf, df)
      
    } ## sim cycle nc
    ## end of simulation
    
    return(adf)
  })
  
  # make graph
  output$graph <- renderPlotly({
    fig <- plot_ly(mySim(), x = ~p, y = ~q, z = ~r, frame = ~t,
                   marker = list(color = ~s, colorscale = "Rainbow", cmin = 0, cmax = 1,size = 15, symbol = "circle", showscale = TRUE))
    fig <- fig %>% add_markers()
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'X'),
                                       yaxis = list(title = 'Y'),
                                       zaxis = list(title = 'Z')))
  }) ## renderPlotly
  
  
} ## server

shinyApp(ui, server)

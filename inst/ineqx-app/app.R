library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(ineqx)
library(ggplot2)
# rsconnect::deployApp("C:/Users/benja/OneDrive - Cornell University/GitHub/ineqx/inst/ineqx-app")

shinyApp(

  ui = dashboardPage(

    title = "Variance decomposition app",
    header = dashboardHeader(
      titleWidth = 350,
      title="Settings"
    ),

    # -------------------------------------------------------------------------------------------- #

    sidebar = dashboardSidebar(
      width = 350,
      collapsed = FALSE,
      materialSwitch("causal", "Causal variance decomposition", status="success"), # toggle
      selectInput("ng", "Number of groups", c(2,3,4), selected=2),
      selectInput("nt", "Number of timepoints", c(1,2,3,4), selected=4),
      uiOutput("tabbox")
    ),

    # -------------------------------------------------------------------------------------------- #

    body = dashboardBody(

      fluidRow(
        column(
          width=12,
          box(
            width=NULL,
            title = "Try the ineqx package!",
            status = "warning",
            solidHeader = TRUE,
            collapsible = TRUE,
            htmlOutput("explanation_input")
          ),
          actionButton("runineqx", "Run ineqx"), # run button
          HTML("<br /><br />"),
          box(
            width=NULL,
            plotOutput("plot_input")
          ),
          htmlOutput("output1"),
          htmlOutput("output2")
        ),
      )
    )
  ),
  server = function(input, output, session) {

    # -------------------------------------------------------------------------------------------- #
    # Input tabbox
    # -------------------------------------------------------------------------------------------- #

    output$tabbox <- renderUI({

      panels <- list()

      for (i in seq_len(input$nt)){


        # Include beta and lambda
        if(isTRUE(input$causal)) {

          tags_mu <- tagList()
          for (j in 1:input$ng) {
            tags_mu[[j]] <- sliderInput(paste0("m",i,j), paste0("Mu ", j), min = -30, max = 30, value = -1+j)
          }
          tags_sigma <- tagList()
          for (j in seq_len(input$ng)) {
            tags_sigma[[j]] <- sliderInput(paste0("s",i,j), paste0("Sigma ", j), min = 6, max = 12, value = 5+j)
          }
          tags_beta <- tagList()
          for (j in seq_len(input$ng)) {
            tags_beta[[j]] <- sliderInput(paste0("b",i,j), paste0("Beta ", j), min = -10, max = 10, value = 1-6*(j==1)+5*(i-1)*(j==1))
          }
          tags_lambda <- tagList()
          for (j in seq_len(input$ng)) {
            tags_lambda[[j]] <- sliderInput(paste0("l",i,j), paste0("Lambda ", j), min = -2, max = 2, value = -3+j+0.5*(i-1)*(j==1), step=0.5)
          }

        } else {
          tags_mu <- tagList()
          for (j in 1:input$ng) {
            tags_mu[[j]] <- sliderInput(paste0("m",i,j), paste0("Mu ", j), min = -30, max = 30, value = 0+10*(j-1)+10*(i-1)*(j==1))
          }
          tags_sigma <- tagList()
          for (j in seq_len(input$ng)) {
            tags_sigma[[j]] <- sliderInput(paste0("s",i,j), paste0("Sigma ", j), min = 6, max = 12, value = 6+1*(j-1)+2*(i-1)*(j==1))
          }
        }

        panels[i] <- tagList(
          tabPanel(
            title=paste0("Time ", i),
            if(isTRUE(input$causal)) {
              fluidRow(
                h4("Pre-treatment inequality", style="color: #000000; text-align: center;"),
                box(
                  title = "Mu",
                  style="color: #000000;",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  collapsed = FALSE,
                  tags_mu
                ),
                box(
                  title = "Sigma",
                  style="color: #000000;",
                  width = 12,
                  status = "success",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  collapsed = FALSE,
                  tags_sigma
                ),
                h4("Treatment effects", style="color: #000000; text-align: center;"),
                box(
                  title = "Beta",
                  style="color: #000000;",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  collapsed = FALSE,
                  tags_beta
                ),
                box(
                  title = "Lambda",
                  style="color: #000000;",
                  width = 12,
                  status = "success",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  collapsed = FALSE,
                  tags_lambda
                )
              )
            } else {
              fluidRow(
                box(
                  title = "Mu",
                  style="color: #000000;",
                  width=12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  collapsed = FALSE,
                  tags_mu
                ),
                box(
                  title = "Sigma",
                  style="color: #000000;",
                  width=12,
                  status = "success",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  collapsed = FALSE,
                  tags_sigma
                )
              )
            }
          )
        )
      }

      do.call(tabBox, args = c(id="tabs", selected="Time 1", width="100%", panels))

    })

    # -------------------------------------------------------------------------------------------- #
    # Input explanation
    # -------------------------------------------------------------------------------------------- #

    output$explanation_input <- renderUI({

      if(isTRUE(input$causal)) {
        HTML("To try the ineqx package, you first choose a few options on the left-hand side and then click on \"Run ineqx\" to run the package.
        <h4>Settings</h4>
        Use the toggle on the left-hand side to decide whether you want to try the descriptive or causal variance decomposition.
        Then choose how many groups and timepoints you want to consider (2 groups and 4 timepoints is a good start).
        <h4>Define inequality across groups and time</h4>
        Next you choose the level of inequality across groups and time.
        Each tab represents one timepoint. For each timepoint and group, you choose the mean and standard deviation of, say, income, which defines the disparities in incomes within and between groups before treatment.
        Then you define the effect of treatment on those means (=beta) and standard deviations (=lambda). While the ineqx package can estimate treatment effects for you, in this example, we assume the true values to be known.
        The plot on the right-hand side visualizes the chosen values.
        As a default, the effect of treatment is 0. If you specify a nonzero treatment effect, you will see that the post-treatment distributions are represented by a dashed distribution.
        <br /><br />
        Then you can go ahead now and run the ineqx package! Note that you have to run the ineqx package again when you change parameters or switch the decomposition type.")
      } else {
        HTML("To try the ineqx package, you first choose a few options on the left-hand side and then click on \"Run ineqx\" to run the package.
        <h4>Settings</h4>
        Use the toggle on the left-hand side to decide whether you want to try the descriptive or causal variance decomposition.
        Then choose how many groups and timepoints you want to consider (2 groups and 4 timepoints is a good start).
        <h4>Define inequality across groups and time</h4>
        Next you choose the level of inequality across groups and time.
        Each tab represents one timepoint. For each timepoint and group, you choose the mean and standard deviation of, say, income, which defines the disparities in incomes within and between groups.
        The plot on the right-hand side visualizes the chosen values.
        <br /><br />
        Then you can go ahead now and run the ineqx package! Note that you have to run the ineqx package again when you change parameters or switch the decomposition type.")
      }

    })

    # -------------------------------------------------------------------------------------------- #
    # Input plot
    # -------------------------------------------------------------------------------------------- #

    output$plot_input <- renderPlot({

      req(input$tabs)
      tab <- as.numeric(gsub("Time ","",input$tabs))
      req(input[[paste0("m",tab,1)]], input[[paste0("s",tab,1)]],
          input[[paste0("m",tab,2)]], input[[paste0("s",tab,2)]])
      if(isTRUE(input$causal)) {
        req(input[[paste0("b",tab,1)]], input[[paste0("l",tab,1)]],
            input[[paste0("b",tab,2)]], input[[paste0("l",tab,2)]])
      }

      # Plot defaults
      p <-
        ggplot(data = data.frame(x = c(-70, 70)), aes(x)) +
        labs(x="", y="", color="", linetype="") +
        scale_x_continuous(breaks = seq(-100,100,10)) +
        scale_y_continuous(limits=c(0, 0.1), breaks = NULL) +
        scale_linetype_manual(values=c("dashed", "solid"), guide = guide_legend(reverse = TRUE)) +
        theme_bw() +
        theme(
          text = element_text(size = 15),
          plot.title = element_text(face="bold", hjust = 0.5),
          axis.title.x = element_text(face="bold"),
          axis.title.y = element_text(face="bold"),
          axis.line = element_line(color = "black"),
          axis.text=element_text(color="black", size = 15),
          panel.grid = element_line(linewidth=0.25),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#A9A9A9", linetype = 2),
          panel.grid.major.x = element_blank(),
          panel.border = element_blank(),
          legend.title = element_text(face="bold"),
          legend.text = element_text(face="bold"),
          legend.position="top",
          legend.margin=margin(0,0,0,0),
          legend.box.margin=margin(0,0,0,0),
          legend.key.width=unit(1.2,"cm"),
          strip.background=element_rect(fill="white", color="white"),
          strip.text=element_text(face="bold", colour = "black", size=rel(1.2)))

      if(isTRUE(input$causal)) {

        p <-
          p +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,1)]], sd = input[[paste0("s",tab,1)]]), linewidth=1, aes(color="Group 1", linetype="(pre)")) +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,1)]] + input[[paste0("b",tab,1)]], sd = input[[paste0("s",tab,1)]] + input[[paste0("l",tab,1)]]), linewidth=1, aes(color="Group 1", linetype="(post)")) +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,2)]], sd = input[[paste0("s",tab,2)]]), linewidth=1, aes(color="Group 2", linetype="(pre)")) +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,2)]] + input[[paste0("b",tab,2)]], sd = input[[paste0("s",tab,2)]] + input[[paste0("l",tab,2)]]), linewidth=1, aes(color="Group 2", linetype="(post)"))

        if(input$ng>2) {
          p <-
            p +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,3)]], sd = input[[paste0("s",tab,3)]]), linewidth=1, aes(color="Group 3", linetype="(pre)")) +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,3)]] + input[[paste0("b",tab,3)]], sd = input[[paste0("s",tab,3)]] + input[[paste0("l",tab,3)]]), linewidth=1, aes(color="Group 3", linetype="(post)"))
        }
        if(input$ng>3) {
          p <-
            p +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,4)]], sd = input[[paste0("s",tab,4)]]), linewidth=1, aes(color="Group 4", linetype="(pre)")) +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,4)]] + input[[paste0("b",tab,4)]], sd = input[[paste0("s",tab,4)]] + input[[paste0("l",tab,4)]]), linewidth=1, aes(color="Group 4", linetype="(post)"))
        }

      } else {

        p <-
          p +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,1)]], sd = input[[paste0("s",tab,1)]]), linewidth=1, aes(color="Group 1")) +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,2)]], sd = input[[paste0("s",tab,2)]]), linewidth=1, aes(color="Group 2"))

        if(input$ng>2) {
          p <-
            p +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,3)]], sd = input[[paste0("s",tab,3)]]), linewidth=1, aes(color="Group 3"))
        }
        if(input$ng>3) {
          p <-
            p +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,4)]], sd = input[[paste0("s",tab,4)]]), linewidth=1, aes(color="Group 4"))
        }

      }

      p

    })

    # -------------------------------------------------------------------------------------------- #
    # Reset outputs when switching between descriptive and causal
    # -------------------------------------------------------------------------------------------- #

    observeEvent(input$causal, {
      output$output1 <- renderUI({ NULL })
      output$output2 <- renderUI({ NULL })
    })

    # -------------------------------------------------------------------------------------------- #
    # Run ineqx (+ all post-run output)
    # -------------------------------------------------------------------------------------------- #

    observeEvent(input$runineqx,{

      set.seed(3)
      ng <- as.integer(input$ng)
      nt <- as.integer(input$nt)

      if(isTRUE(input$causal)) {

        # Build parameter data.frame from slider inputs
        rows <- list()
        for(i in seq_len(nt)) {
          for(j in seq_len(ng)) {
            rows[[length(rows) + 1]] <- data.frame(
              group  = j,
              time   = i,
              pi     = 1 / ng,
              mu0    = input[[paste0("m",i,j)]],
              sigma0 = input[[paste0("s",i,j)]],
              beta   = input[[paste0("b",i,j)]],
              lambda = input[[paste0("l",i,j)]] / input[[paste0("s",i,j)]]  # convert additive to log-scale
            )
          }
        }
        param_df <- do.call(rbind, rows)

        params <- ineqx_params(data = param_df, ref = 1)

        if (params$type == "cross_sectional") {
          ineqx.out <- ineqx(params = params, se = "none")
        } else {
          ineqx.out <- ineqx(params = params, order = "shapley", se = "none")
        }

      } else {

        # Build data for descriptive decomposition
        n <- 1000
        dat_list <- list()
        for(i in seq_len(nt)) {
          for(j in seq_len(ng)) {
            dat_list[[length(dat_list) + 1]] <- data.frame(
              id = seq_len(n),
              time = i,
              group = j,
              y = rnorm(n, input[[paste0("m",i,j)]], input[[paste0("s",i,j)]])
            )
          }
        }
        dat <- do.call(rbind, dat_list)

        ineqx.out <- ineq(
          y = "y", group = "group", time = "time",
          data = dat, ref = 1, ystat = "Var"
        )
      }

      # Output 1: W/B levels over time --------------------------------------------------------- #

      output$output1 <- renderUI({
        box(
          width=NULL,
          tagList(
            HTML("<h4>Within- and between-group inequality over time</h4>
               This plot displays the development of within-group, between-group, and total inequality over time."),
            plotOutput("plot_wibe"),
          )
        )
      })

      output$plot_wibe <- renderPlot({
        if(isTRUE(input$causal)) {
          # Causal decomposition plot (bar chart for cross-sectional, line for longit/shapley)
          plot(ineqx.out)
        } else {
          plot(ineqx.out, type = "wibe")
        }
      })

      # Output 2: Decomposition details -------------------------------------------------------- #

      output$output2 <- renderUI({
        if(!isTRUE(input$causal)) {
          box(
            width=NULL,
            tagList(
              HTML("<h4>Decomposition of change over time</h4>
               This plot decomposes the change in inequality relative to the reference period:
               <ul>
                <li><b>Means (mu)</b>: How much did inequality change due to changes in group means?</li>
                <li><b>Dispersions (sigma)</b>: How much did inequality change due to changes in within-group dispersions?</li>
                <li><b>Composition (pi)</b>: How much did inequality change due to changes in group composition?</li>
               </ul>"),
              plotOutput("plot_ineqx")
            )
          )
        }
      })

      output$plot_ineqx <- renderPlot({
        plot(ineqx.out, type = "deltas")
      })

    })

  }
)

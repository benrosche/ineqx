library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(DT)
library(tidyr)
library(ineqx)
library(ggplot2)
library(ggthemes)

shinyApp(

  ui = dashboardPage(

    title = "Variance decomposition app",
    header = dashboardHeader(title="Settings"),

    # -------------------------------------------------------------------------------------------- #

    sidebar = dashboardSidebar(
      collapsed = FALSE,

      materialSwitch("causal", "Causal variance decomposition", status="success"),
      selectInput("ng", "Number of groups", c(2,3,4), selected=2),
      selectInput("nt", "Number of timepoints", c(1,2,3,4), selected=3)
    ),

    # -------------------------------------------------------------------------------------------- #

    body = dashboardBody(

      fluidRow(
        column(
          width=7,
          box(
            width=NULL,
            plotOutput("dist")
          ),
          box(
            width=NULL,
            shinyjs::useShinyjs(),  # Set up shinyjs
            disabled(actionButton("runineqx", "Run ineqx")),
            htmlOutput("explanation_output1"),
            plotOutput("ineqxplot")
          )
        ),
        column(
          width=5,
          uiOutput("tabbox")
        )

      ),
      fluidRow(
        box(
          width=12,
          htmlOutput("explanation_output2"),
          DT::dataTableOutput("ineqxtable")
        )
      )
    )
  ),
  server = function(input, output, session) {

    # -------------------------------------------------------------------------------------------- #
    # Tabbox for input on the right-hand side
    # -------------------------------------------------------------------------------------------- #

    output$tabbox <- renderUI({

      panels <- list()

      panels[1] <- tagList(
        tabPanel(
          title="Time",
          htmlOutput("explanation_input")
        )
      )

      for (i in seq_len(input$nt)){

        tags_mu <- tagList()
        for (j in 1:input$ng) {
          tags_mu[[j]] <- sliderInput(paste0("m",i,j), paste0("Mu ", j), min = -50, max = 50, value = 0+10*j)
        }
        tags_sigma <- tagList()
        for (j in seq_len(input$ng)) {
          tags_sigma[[j]] <- sliderInput(paste0("s",i,j), paste0("Sigma ", j), min = 5, max = 10, value = 0+1*j)
        }
        # Include beta and lambda
        if(isTRUE(input$causal)) {
          tags_beta <- tagList()
          for (j in seq_len(input$ng)) {
            tags_beta[[j]] <- sliderInput(paste0("b",i,j), paste0("Beta ", j), min = -10, max = 10, value = 0)
          }
          tags_lambda <- tagList()
          for (j in seq_len(input$ng)) {
            tags_lambda[[j]] <- sliderInput(paste0("l",i,j), paste0("Lambda ", j), min = -1, max = 1, value = 0)
          }
        }

        panels[i+1] <- tagList(
          tabPanel(
            title=i,
            if(isTRUE(input$causal)) {
              fluidRow(
                h4("Pre-treatment inequality"),
                box(
                  title = "Mu",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  tags_mu
                ),
                box(
                  title = "Sigma",
                  width = 12,
                  status = "success",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  tags_sigma
                ),
                h4("Treatment effect"),
                box(
                  title = "Beta",
                  width = 12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  tags_beta
                ),
                box(
                  title = "Lambda",
                  width = 12,
                  status = "success",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  tags_lambda
                )
              )
            } else {
              fluidRow(
                box(
                  title = "Mu",
                  width=12,
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  tags_mu
                ),
                box(
                  title = "Sigma",
                  width=12,
                  status = "success",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  tags_sigma
                )
              )
            }
          )
        )
      }

      do.call(tabBox, args = c(id="tabs", selected="Time", width="100%", panels))

    })


    # -------------------------------------------------------------------------------------------- #
    # Input plot in the middle
    # -------------------------------------------------------------------------------------------- #

    output$dist <- renderPlot({

      if(isTRUE(input$causal)) {

        p <- ggplot(data = data.frame(x = c(-90, 90)), aes(x)) +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",input$tabs,1)]], sd = input[[paste0("s",input$tabs,1)]]), size=1, color="red") +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",input$tabs,1)]] + input[[paste0("b",input$tabs,1)]], sd = input[[paste0("s",input$tabs,1)]] + input[[paste0("l",input$tabs,1)]]), size=1, color="red", linetype="dashed") +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",input$tabs,2)]], sd = input[[paste0("s",input$tabs,2)]]), size=1, color="blue") +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",input$tabs,2)]] + input[[paste0("b",input$tabs,2)]], sd = input[[paste0("s",input$tabs,2)]] + input[[paste0("l",input$tabs,2)]]), size=1, color="blue", linetype="dashed") +
          labs(x="", y="") + scale_x_continuous(breaks = seq(-90,90,15)) + scale_y_continuous(limits=c(0, 0.1), breaks = NULL) + theme_economist()

        if(input$ng>2) {
          p <-
            p +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",input$tabs,3)]], sd = input[[paste0("s",input$tabs,3)]]), size=1, color="green") +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",input$tabs,3)]] + input[[paste0("b",input$tabs,3)]], sd = input[[paste0("s",input$tabs,3)]] + input[[paste0("l",input$tabs,3)]]), size=1, color="green" , linetype="dashed")
        }
        if(input$ng>3) {
          p <-
            p +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",input$tabs,4)]], sd = input[[paste0("s",input$tabs,4)]]), size=1, color="yellow") +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",input$tabs,4)]] + input[[paste0("b",input$tabs,4)]], sd = input[[paste0("s",input$tabs,4)]] + input[[paste0("l",input$tabs,4)]]), size=1, color="yellow", linetype="dashed")
        }

      } else {

        p <- ggplot(data = data.frame(x = c(-90, 90)), aes(x)) +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",input$tabs,1)]], sd = input[[paste0("s",input$tabs,1)]]), size=1, color="red") +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",input$tabs,2)]], sd = input[[paste0("s",input$tabs,2)]]), size=1, color="blue") +
          labs(x="", y="") + scale_x_continuous(breaks = seq(-90,90,15)) + scale_y_continuous(limits=c(0, 0.1), breaks = NULL) + theme_economist()

        if(input$ng>2) {
          p <- p + stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",input$tabs,3)]], sd = input[[paste0("s",input$tabs,3)]]), size=1, color="green")
        }
        if(input$ng>3) {
          p <- p + stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",input$tabs,4)]], sd = input[[paste0("s",input$tabs,4)]]), size=1, color="yellow")
        }

      }

      p

    })

    # -------------------------------------------------------------------------------------------- #
    # Input explanation
    # -------------------------------------------------------------------------------------------- #

    output$explanation_input <- renderUI({
      if(isTRUE(input$causal)) {
        HTML("<h4>Try the ineqx package!</h4>
        To try the ineqx package, you have to choose a few options on the left and right-hand side and you can click on \"Run ineqx\" to run the package.
        <h4>Settings</h4>
        Use the toggle on the left-hand side to decide whether you want to try the descriptive or causal variance decomposition.
        Also choose how many groups and timepoints you want to consider. 3 groups and 4 timepoints is a good start.
        <h4>Define inequality across groups and time</h4>
        Next you have to choose inequality across groups and time. You this is in the other tabs here on the right-hand side.
        In tab 1 you define the inequality within and across groups at timepoint 1, in tab 2 you define inequality at timepoint 2, and so on.
        <h4>Causal variance decomposition</h4>
        To try the causal variance decomposition, you define the pre-treatment mean and standard deviation of, say, income for each group at each timepoint.
        Then you also define the effect of treatment on those means (=beta) and standard deviations (=lambda).
        The plot in the middle visualizes the chosen values. As a default, the treatment effect is 0.
        Once you specify some values instead, you will see that inequality post treatment is visualized by a dashed distribution.
        <br /><br />
        <b>Okay, now you can go ahead and run the ineqx package!</b>")
      } else {
        HTML("<h4>Try the ineqx package!</h4>
        To try the ineqx package, you have to choose a few options on the left and right-hand side and you can click on \"Run ineqx\" to run the package.
        <h4>Settings</h4>
        Use the toggle on the left-hand side to decide whether you want to try the descriptive or causal variance decomposition.
        Also choose how many groups and timepoints you want to consider. 3 groups and 4 timepoints is a good start.
        <h4>Define inequality across groups and time</h4>
        Next you have to choose inequality across groups and time. You this is in the other tabs here on the right-hand side.
        In tab 1 you define the inequality within and across groups at timepoint 1, in tab 2 you define inequality at timepoint 2, and so on.
        <h4>Descriptive variance decomposition</h4>
        To try the descriptive variance decomposition, you just define the mean and standard deviation of, say, income of each group at each timepoint.
        The plot in the middle visualizes the chosen values. I have prespecified some values for you that are a good starting point.
        <br /><br />
        <b>Okay, now you can go ahead and run the ineqx package!</b>")
      }

    })

    # -------------------------------------------------------------------------------------------- #
    # Enable action button
    # -------------------------------------------------------------------------------------------- #

    observeEvent(input$tabs, {
      if(!input$tabs=="Time") {
        shinyjs::enable('runineqx')
      }
    })

    # -------------------------------------------------------------------------------------------- #
    # Run ineqx
    # -------------------------------------------------------------------------------------------- #

    observeEvent(input$runineqx,{

      set.seed(1)
      n <- 100

      dat <-
        tibble() %>%
        tidyr:::expand(id=seq_len(n), time=seq_len(input$nt), group=seq_len(input$ng), x=c(0,1), y=NA)

      if(isTRUE(input$causal)) {
        for(i in seq_len(input$nt)) {
          for(j in seq_len(input$ng)) {
            dat[dat$time==i & dat$group==j & dat$x==0,] <- dat %>% dplyr::filter(time==i & group==j & x==0) %>% dplyr::mutate(y=rnorm(n, input[[paste0("m",i,j)]], input[[paste0("s",i,j)]]))
            dat[dat$time==i & dat$group==j & dat$x==1,] <- dat %>% dplyr::filter(time==i & group==j & x==1) %>% dplyr::mutate(y=rnorm(n, input[[paste0("m",i,j)]]+input[[paste0("b",i,j)]], input[[paste0("s",i,j)]]+input[[paste0("l",i,j)]]))
          }
        }
      } else {
        for(i in seq_len(input$nt)) {
          for(j in seq_len(input$ng)) {
            dat[dat$time==i & dat$group==j,] <- dat %>% dplyr::filter(time==i & group==j) %>% dplyr::mutate(y=rnorm(2*n, input[[paste0("m",i,j)]], input[[paste0("s",i,j)]]))
          }
        }
      }

      if(isTRUE(input$causal)) {

        AME_mu <- list()
        AME_mu[[1]] <- tibble(time=numeric(), group=numeric(), effect=numeric())
        AME_mu[[2]] <- tibble(time=numeric(), effect=numeric())

        AME_sigma <- list()
        AME_sigma[[1]] <- tibble(time=numeric(), group=numeric(), effect=numeric())
        AME_sigma[[2]] <- tibble(time=numeric(), effect=numeric())

        for(i in seq_len(input$nt)) {
          for(j in seq_len(input$ng)) {
            AME_mu[[1]] <- AME_mu[[1]] %>% add_row(time=i, group=j, effect=input[[paste0("b",i,j)]])
            AME_sigma[[1]] <- AME_sigma[[1]] %>% add_row(time=i, group=j, effect=input[[paste0("l",i,j)]])
          }
          AME_mu[[2]] <- AME_mu[[1]] %>% group_by(time) %>% dplyr::summarise(time=mean(time), effect=mean(effect))
          AME_sigma[[2]] <- AME_sigma[[1]] %>% group_by(time) %>% dplyr::summarise(time=mean(time), effect=mean(effect))
        }

      }

      if(isTRUE(input$causal)) {
        ineqx.out <- ineqx(treat="x", y="y", group="group", time="i.time", AME_mu=AME_mu, AME_sigma=AME_sigma, ref=1, dat=dat)
      } else {
        ineqx.out <- ineqx(y="y", group="group", time="i.time", ref=1, dat=dat)
      }

      output$ineqxplot <- renderPlot({ plot(ineqx.out, type="dT") })

      output$ineqxtable <- DT::renderDataTable(DT::datatable(ineqx.out$dT[[1]] %>% dplyr::mutate(across(where(is.numeric), round, digits=2)), options = list(pageLength = 5, dom = 'tip'), selection = "none", rownames = FALSE))

      output$explanation_output1 <- renderUI({
        if(isTRUE(input$causal)) {
          HTML("<h4>Interpreting the results</h4>
               This plot displays the results of the causal variance decomposition.
               The change of four quantities over time are depicted:
               <ul>
                <li>The between-group effect ...</li>
                <li>The within-group effect ...</li>
                <li>The compositional-group effect ...</li>
                <li>The pre-treatment effect ...</li>
               </ul>
               Below you can find the output in table form.
               <br /><br />")
        } else {
          HTML("<h4>Interpreting the results</h4>
               This plot displays the results of the descriptive variance decomposition.
               The change of three quantities over time are depicted:
               <ul>
                <li>The between-group effect ...</li>
                <li>The within-group effect ...</li>
                <li>The compositional-group effect ...</li>
               </ul>
               Below you can find the output in table form.
               <br /><br />")
        }
        })

      output$explanation_output2 <- renderUI({
        if(isTRUE(input$causal)) {
          HTML("<h4>Interpreting the results</h4>
               This table displays the results of the causal variance decomposition.
               The change of four quantities over time are depicted:
               <ul>
                <li>The between-group effect ...</li>
                <li>The within-group effect ...</li>
                <li>The compositional-group effect ...</li>
                <li>The pre-treatment effect ...</li>
               </ul>")
        } else {
          HTML("<h4>Interpreting the results</h4>
               This table displays the results of the causal variance decomposition.
               The change of three quantities over time are depicted:
               <ul>
                <li>The between-group effect ...</li>
                <li>The within-group effect ...</li>
                <li>The compositional-group effect ...</li>
                <li>The pre-treatment effect ...</li>
               </ul>")
        }
      })

    })

  }

)


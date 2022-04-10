library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(DT)
library(tidyr)
library(ineqx)
library(ggplot2)
library(ggthemes)
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
          htmlOutput("output2"),
          htmlOutput("output3")
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

      tab <- as.numeric(gsub("Time ","",input$tabs))

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
          panel.grid = element_line(size=0.25),
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
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,1)]], sd = input[[paste0("s",tab,1)]]), size=1, aes(color="Group 1", linetype="(pre)")) +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,1)]] + input[[paste0("b",tab,1)]], sd = input[[paste0("s",tab,1)]] + input[[paste0("l",tab,1)]]), size=1, aes(color="Group 1", linetype="(post)")) +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,2)]], sd = input[[paste0("s",tab,2)]]), size=1, aes(color="Group 2", linetype="(pre)")) +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,2)]] + input[[paste0("b",tab,2)]], sd = input[[paste0("s",tab,2)]] + input[[paste0("l",tab,2)]]), size=1, aes(color="Group 2", linetype="(post)"))

        if(input$ng>2) {
          p <-
            p +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,3)]], sd = input[[paste0("s",tab,3)]]), size=1, aes(color="Group 3", linetype="(pre)")) +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,3)]] + input[[paste0("b",tab,3)]], sd = input[[paste0("s",tab,3)]] + input[[paste0("l",tab,3)]]), size=1, aes(color="Group 3", linetype="(post)"))
        }
        if(input$ng>3) {
          p <-
            p +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,4)]], sd = input[[paste0("s",tab,4)]]), size=1, aes(color="Group 4", linetype="(pre)")) +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,4)]] + input[[paste0("b",tab,4)]], sd = input[[paste0("s",tab,4)]] + input[[paste0("l",tab,4)]]), size=1, aes(color="Group 4", linetype="(post)"))
        }

      } else {

        p <-
          p +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,1)]], sd = input[[paste0("s",tab,1)]]), size=1, aes(color="Group 1")) +
          stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,2)]], sd = input[[paste0("s",tab,2)]]), size=1, aes(color="Group 2"))

        if(input$ng>2) {
          p <-
            p +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,3)]], sd = input[[paste0("s",tab,3)]]), size=1, aes(color="Group 3"))
        }
        if(input$ng>3) {
          p <-
            p +
            stat_function(fun = dnorm, n = 10000, args = list(mean = input[[paste0("m",tab,4)]], sd = input[[paste0("s",tab,4)]]), size=1, aes(color="Group 4"))
        }

      }

      p

    })

    # -------------------------------------------------------------------------------------------- #
    # Run ineqx (+ all post-run output)
    # -------------------------------------------------------------------------------------------- #

    observeEvent(input$runineqx,{

      set.seed(3)

      if(isTRUE(input$causal)) {

        n <- 10000

        dat <-
          tibble() %>%
          tidyr:::expand(id=seq_len(n), time=seq_len(input$nt), group=seq_len(input$ng), x=c(0,1), t=c(0,1), y=NA)

        for(i in seq_len(input$nt)) {
          for(j in seq_len(input$ng)) {
            dat[dat$time==i & dat$group==j & dat$t==0,] <- dat %>% dplyr::filter(time==i & group==j & t==0) %>% dplyr::mutate(y=rnorm(2*n, input[[paste0("m",i,j)]], input[[paste0("s",i,j)]]))
            dat[dat$time==i & dat$group==j & dat$t==1 & dat$x==1,] <- dat %>% dplyr::filter(time==i & group==j & x==1 & t==1) %>% dplyr::mutate(y=rnorm(n, input[[paste0("m",i,j)]]+input[[paste0("b",i,j)]], input[[paste0("s",i,j)]]+input[[paste0("l",i,j)]]))
          }
        }
      } else {

        n <- 1000

        dat <-
          tibble() %>%
          tidyr:::expand(id=seq_len(n), time=seq_len(input$nt), group=seq_len(input$ng), y=NA)

        for(i in seq_len(input$nt)) {
          for(j in seq_len(input$ng)) {
            dat[dat$time==i & dat$group==j,] <- dat %>% dplyr::filter(time==i & group==j) %>% dplyr::mutate(y=rnorm(n, input[[paste0("m",i,j)]], input[[paste0("s",i,j)]]))
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

        ineqx.out <- ineqx(treat="x", post="t", y="y", ystat="Var", group="group", time="i.time", AME_mu=AME_mu, AME_sigma=AME_sigma, ref=1, dat=dat)

        wibe.out  <-
          wibe(y="y", group="group", time="i.time", long=T, dat=
                 dat %>%
                 dplyr::filter(x==0))[[2]] %>%
          dplyr::filter(variable %in% c("VarW", "VarB", "VarT")) %>%
          dplyr::mutate(x=0) %>%
          add_row(
            wibe(y="y", group="group", time="i.time", long=T, dat=
                   dat %>%
                   dplyr::filter(x==1))[[2]] %>%
              dplyr::filter(variable %in% c("VarW", "VarB", "VarT")) %>%
              dplyr::mutate(x=1)
          )

      } else {

        ineqx.out <- ineqx(y="y", ystat="Var", group="group", time="i.time", ref=1, dat=dat)

        wibe.out  <-
          wibe(y="y", group="group", time="i.time", long=T, dat=dat)[[2]] %>%
          dplyr::filter(variable %in% c("VarW", "VarB", "VarT"))

      }

      # Output 1 --------------------------------------------------------------------------------- #

      output$output1 <- renderUI({

        box(
          width=NULL,
          tagList(
            HTML("<h4>Within- and between-group inequality over time</h4>
               This plot displays the development of within-group, between-group, and total inequality over time in absolute values."),
            plotOutput("plot_wibe"),
          )
        )

      })

      output$plot_wibe <- renderPlot({

        if(isTRUE(input$causal)) {
          p <- ggplot(aes(x=time, y=value, color=factor(variable, levels = c("VarB", "VarW", "VarT")), linetype=factor(x, levels = c(0,1), labels = c("pre", "post"))), data = wibe.out)
        } else {
          p <- ggplot(aes(x=time, y=value, color=factor(variable, levels = c("VarB", "VarW", "VarT"))), data = wibe.out)
        }

        p <-
          p +
          geom_line() +
          labs(x="Time", y="Variance", color="", linetype="") +
          scale_color_manual(values=c("#F8766D", "#00BFC4", "#000000")) +
          scale_linetype_manual(values = c("solid", "dashed")) +
          scale_x_continuous(breaks=seq(1:4)) +
          #scale_y_continuous(limits = c(0, 500)) +
          theme_bw() +
          theme(
            text = element_text(size = 15),
            plot.title = element_text(face="bold", hjust = 0.5),
            axis.title.x = element_text(face="bold"),
            axis.title.y = element_text(face="bold"),
            axis.line = element_line(color = "black"),
            axis.text=element_text(color="black", size = 15),
            panel.grid = element_line(size=0.25),
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

        p

      })

      # Output 2 --------------------------------------------------------------------------------- #

      output$output2 <- renderUI({

        box(
          width=NULL,
          tagList(
            if(isTRUE(input$causal)) {
              HTML("<h4>Decomposition results</h4>
               This plot displays the results of the causal variance decomposition.
               The change of four quantities are visualized:
               <ul>
                <li><b>The between-group effect</b>: How much did the treatment effect on the variance change due to changes in the effect of treatment on between-group inequality?</li>
                <li><b>The within-group effect</b>: How much did the treatment effect on the variance change due to changes in the effect of treatment on within-group inequality?</li>
                <li><b>The compositional-group effect</b>: How much did the treatment effect on the variance change due to changes in the composition of groups? This value is zero because groups sizes are constant in the example.</li>
                <li><b>The pre-treatment effect</b>: How much did the treatment effect on the variance change due to changes in pre-treatment inequality? (Note that this quantity can fluctuate due to the probabilistic nature of the distributions.)</li>
               </ul>")
            } else {
              HTML("<h4>Decomposition results</h4>
               This plot displays the results of the descriptive variance decomposition.
               The change of three quantities are visualized:
               <ul>
                <li><b>The between-group effect</b>: How much did the total variance change due to changes in between-group inequality?</li>
                <li><b>The within-group effect</b>: How much did the total variance change due to changes in within-group inequality?</li>
                <li><b>The compositional-group effect</b>: How much did the total variance change due to changes in the composition of groups? This value is zero because groups sizes are constant in the example.</li>
               </ul>")
            },
            plotOutput("plot_ineqx")
          )
        )

      })

      output$plot_ineqx <- renderPlot({ plot(ineqx.out, type="dT") + scale_x_continuous(breaks=seq(1,4,1)) })

      # Output 3 --------------------------------------------------------------------------------- #

      output$output3 <- renderUI({

        box(
          width=NULL,
          tagList(
            HTML("<h4>The ineqx package offers different table and plot types to facilitate interpretation</h4>"),
            htmlOutput("select"),
            selectInput(
              inputId="select_output3",
              label="Choose from list to see other output types",
              choices = list(
                "Table \"Total\"" = 1,
                "Plot \"Total\"" = 2,
                "Plot \"Shares\"" = 3,
                "Plot \"Treatment effect on Mu\"" = 4,
                "Plot \"Treatment effect on Sigma\"" = 5),
              selected = 1)
          )
        )
      })

      observeEvent(
        input$select_output3, {
          output$select <- renderUI({
            if(input$select_output3==1) {
              tagList(
                HTML("This table displays the same output as the plot above."),
                DT::renderDataTable(DT::datatable(ineqx.out$dT[[1]] %>% dplyr::mutate(across(where(is.numeric), round, digits=2)), options = list(pageLength = 5, dom = 'tip'), selection = "none", rownames = FALSE))
              )
            } else if(input$select_output3==2){
              tagList(
                HTML("This plot displays the same output as the plot above as well as the change in the variance, dVarT (yellow).
                     We can see that the total change in the four components (dT) equals the change in the variance (dVarT). (The overlay is not perfect due to the probabilistic nature of the data.)
                     The reason why dT and dVarT are equal is that the true treatment effect values are known and no other forces influence the variance in this example.
                     In practice, the more the changes in the variance are driven by changes in the treatment effect, the closer the two lines will follow each other."),
                renderPlot({ plot(ineqx.out, type="dPA") + scale_x_continuous(breaks=seq(1,4,1)) })
              )
            } else if(input$select_output3==3){
              tagList(
                HTML("This plot displays the size of within-group, between-group, compositional, and pre-treatment effects are relative to each other in absolute values."),
                renderPlot({ plot(ineqx.out, type="dTS") })
              )
            } else if(input$select_output3==4){
              tagList(
                HTML("This plot displays the treatment effects on the mean. In a real application those values would be estimated."),
                renderPlot({ plot(ineqx.out, type="dMuP") + labs(color="", group="", linetype="") })
              )
            } else if(input$select_output3==5){
              tagList(
                HTML("This plot displays the treatment effects on the variance. In a real application those values would be estimated."),
                renderPlot({ plot(ineqx.out, type="dSigmaP") + labs(color="", group="", linetype="") })
              )
            }
          })
        })

    })

  }
)


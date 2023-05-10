using Genie, Genie.Renderer, Genie.Renderer.Html, Genie.Renderer.Json, Genie.Requests, SimpleWebsockets

include("src/solve.jl")

Genie.config.run_as_server = true
Genie.config.cors_headers["Access-Control-Allow-Origin"] = "https://main.d159dyqns1i521.amplifyapp.com/"
# This has to be this way - you should not include ".../*"
Genie.config.cors_headers["Access-Control-Allow-Headers"] = "Content-Type"
Genie.config.cors_headers["Access-Control-Allow-Methods"] ="GET,POST,PUT,DELETE,OPTIONS" 
Genie.config.cors_allowed_origins = ["*"]

route("/") do
  html("Hello World")
end

server = WebsocketServer()

@async serve(server, 8080, "teemaserver.cloud")

listen(server, :client) do client 
  route("/solving" ,method="POST") do 
      return JSON.json(doSolving(jsonpayload()["mesherOutput"], jsonpayload()["solverInput"], jsonpayload()["solverAlgoParams"], client))
  end
end
using HTTP

# insert your real key here!
access_key = "your_key"

headers = Dict("X-Yandex-Weather-Key"=>access_key)

query = "{
  weatherByPoint(request: { lat: 52.37125, lon: 4.89388 }) {
    now {
      temperature
    }
  }
}"
data = Dict("query"=> query)
res = HTTP.request("POST","https://api.weather.yandex.ru/graphql/query", headers, data)
res = JSON.parse(String(res.body));
print(response.content)
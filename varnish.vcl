vcl 4.0;

acl purgers {
     "127.0.0.1";
}

backend default {
    .host = "127.0.0.1";
    .port = "8080";
}

sub vcl_recv {
    if (req.method == "PURGE") {
      if (!client.ip ~ purgers) {
        return (synth(405));
      }
      return (purge);
    }

    unset req.http.Cookie;
}

sub vcl_backend_response {
    if (beresp.status == 404 || beresp.status == 301 || beresp.status == 500) {
      set beresp.ttl = 15s;
      return (deliver);
    }

    set beresp.ttl = 1d;
    set beresp.grace = 30s;

    unset beresp.http.Server;
    unset beresp.http.Vary;
}

sub vcl_deliver {
    unset resp.http.Server;
    unset resp.http.Via;
    unset resp.http.X-Varnish;
}
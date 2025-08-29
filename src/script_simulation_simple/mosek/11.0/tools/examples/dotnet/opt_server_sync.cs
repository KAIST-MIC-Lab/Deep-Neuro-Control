/*
  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.

  File:    opt_server_sync.cs

  Purpose :   Demonstrates how to use MOSEK OptServer
              to solve optimization problem synchronously
*/
using System;

namespace mosek.example
{
  class msgclass : mosek.Stream
  {
    public override void streamCB (string msg)
    {
      Console.Write ("{0}", msg);
    }
  }

  public class simple
  {
    public static void Main (string[] args)
    {
      if (args.Length == 0)
      {
        Console.WriteLine ("Missing arguments, syntax is:");
        Console.WriteLine ("  opt_server_sync inputfile addr [certpath]");
      }
      else
      {
        String inputfile = args[0];
        String addr      = args[1];
        String cert      = args.Length < 3 ? null : args[2];

        using (mosek.Task task = new mosek.Task())
        {
          task.set_Stream (mosek.streamtype.log, new msgclass ());

          // Load some data into the task
          task.readdata (inputfile);

          // Set OptServer URL
          task.putoptserverhost(addr);

          // Path to certificate, if any
          if (cert != null)
            task.putstrparam(sparam.remote_tls_cert_path, cert);

          // Optimize remotely, no access token
          mosek.rescode trm = task.optimize ();

          task.solutionsummary (mosek.streamtype.log);
        }
      }
    }
  }
}

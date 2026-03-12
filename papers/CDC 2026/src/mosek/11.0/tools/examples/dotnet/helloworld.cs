////
//  Copyright: Copyright (c) MOSEK ApS, Denmark. All rights reserved.
//
//  File:      helloworld.cs
//
//  The most basic example of how to get started with MOSEK.

using mosek;
using System;

public class helloworld {
  public static void Main() {

    using (Task task = new Task()) {        // Create Task

      task.appendvars(1);                          // 1 variable x
      task.putcj(0, 1.0);                          // c_0 = 1.0
      task.putvarbound(0, boundkey.ra, 2.0, 3.0);  // 2.0 <= x <= 3.0
      task.putobjsense(objsense.minimize);         // minimize

      task.optimize();                      // Optimize

      double[] x = task.getxx(soltype.itr);        // Get solution
      Console.WriteLine("Solution x = " + x[0]);   // Print solution
    }
  }
}

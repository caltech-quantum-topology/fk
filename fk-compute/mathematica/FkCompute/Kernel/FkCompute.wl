BeginPackage["FkCompute`"]

FkCompute::usage =
  "FkCompute[braid, degree, opts] computes the FK invariant using the installed fkcompute Python package.\n" <>
  "FkCompute[configPath] runs from a YAML/JSON config file (same format as the fk CLI).";

FkComputeVersion::usage =
  "FkComputeVersion[] returns version info for Python and fkcompute as seen by the wrapper.";

$FkComputePythonExecutable::usage =
  "$FkComputePythonExecutable sets the default Python executable used by FkCompute (e.g. \"python3\" or a full path).";

SetFkComputePythonExecutable::usage =
  "SetFkComputePythonExecutable[exe_String] sets the Python executable used by FkCompute.";

GetFkComputePythonExecutable::usage =
  "GetFkComputePythonExecutable[] returns the current Python executable used by FkCompute.";

FkCompute::pyfail = "Python bridge failed (exit code `1`). stderr: `2`";
FkCompute::pybadjson = "Python bridge returned non-JSON output: `1`";
FkCompute::pyerror = "fkcompute raised `1`: `2`";

Begin["`Private`"]

$FkComputePythonExecutable = Automatic;

SetFkComputePythonExecutable[exe_String] := ($FkComputePythonExecutable = exe)

GetFkComputePythonExecutable[] := $FkComputePythonExecutable

Options[FkCompute] = {
  "PythonExecutable" -> Automatic,
  "Symbolic" -> False,
  "Threads" -> Automatic,
  "Verbose" -> False,
  "Name" -> Automatic,
  "MaxWorkers" -> Automatic,
  "ChunkSize" -> Automatic,
  "IncludeFlip" -> Automatic,
  "MaxShifts" -> Automatic,
  "SaveData" -> Automatic,
  "SaveDir" -> Automatic,
  "ILPFile" -> Automatic,
  "InversionFile" -> Automatic,
  "PartialSigns" -> Automatic,
  "Preset" -> Automatic
};

Options[FkComputeVersion] = {"PythonExecutable" -> Automatic};

FkComputePythonExe[opt_] := Replace[
  opt,
  Automatic :> Replace[$FkComputePythonExecutable, Automatic :> "python3"]
];

FkComputeDropAutomatic[assoc_Association] := KeySelect[assoc, assoc[#] =!= Automatic &];

FkComputeToKwargs[opts_Association] := Module[{m},
  m = <|
    "symbolic" -> opts["Symbolic"],
    "threads" -> opts["Threads"],
    "verbose" -> opts["Verbose"],
    "name" -> opts["Name"],
    "max_workers" -> opts["MaxWorkers"],
    "chunk_size" -> opts["ChunkSize"],
    "include_flip" -> opts["IncludeFlip"],
    "max_shifts" -> opts["MaxShifts"],
    "save_data" -> opts["SaveData"],
    "save_dir" -> opts["SaveDir"],
    "ilp_file" -> opts["ILPFile"],
    "inversion_file" -> opts["InversionFile"],
    "partial_signs" -> opts["PartialSigns"],
    "preset" -> opts["Preset"]
  |>;
  FkComputeDropAutomatic[m]
];

FkComputeCallBridge[pyExe_, req_Association] := Module[
  {reqJson, proc, code, out, errout, parsed, err, mode},

  mode = Lookup[req, "mode", Automatic];
  reqJson = ExportString[req, "RawJSON"];

  (* Run bridge directly (no shell). Use the real python binary path, e.g. .../bin/python3.13 *)
  proc = RunProcess[{pyExe, "-m", "fkcompute.mathematica_bridge"}, All, reqJson];

  If[proc === $Failed,
    Message[FkCompute::pyfail, $Failed, "RunProcess returned $Failed (python executable not runnable)."];
    Return[$Failed];
  ];

  code = proc["ExitCode"];
  out = proc["StandardOutput"];
  errout = proc["StandardError"];

  parsed = Quiet@Check[ImportString[out, "RawJSON"], $Failed];
  If[parsed === $Failed,
    Message[FkCompute::pybadjson, StringTake[ToString[out, InputForm], UpTo[800]]];
    Message[FkCompute::pyfail, code, StringTake[ToString[errout, InputForm], UpTo[800]]];
    Return[$Failed];
  ];

  If[TrueQ[parsed["ok"]],
    (* Convert compute results from string -> WL expression; leave other modes unchanged *)
    If[mode === "compute",
      Replace[parsed["result"],
        s_String :> Module[{t},
          t = StringReplace[s, {"\\\r\n" -> "", "\\\n" -> "", "\r" -> " ", "\n" -> " "}];
          ToExpression[t, InputForm]
        ]
      ],
      parsed["result"]
    ],
    err = parsed["error"];
    Message[FkCompute::pyerror, Lookup[err, "type", "Error"], Lookup[err, "message", ""]];
    If[code =!= 0 && StringLength[ToString[errout]] > 0,
      Message[FkCompute::pyfail, code, StringTake[ToString[errout, InputForm], UpTo[800]]]
    ];
    $Failed
  ]
];

FkCompute[braid_List, degree_Integer, opts : OptionsPattern[]] := Module[
  {pyExe, optsAssoc, req},
  pyExe = FkComputePythonExe[OptionValue["PythonExecutable"]];
  optsAssoc = Join[
    Association[Options[FkCompute]],
    Association@FilterRules[{opts}, Options[FkCompute]]
  ];
  req = <|
    "mode" -> "compute",
    "braid" -> braid,
    "degree" -> degree,
    "kwargs" -> FkComputeToKwargs[optsAssoc]
  |>;
  FkComputeCallBridge[pyExe, req]
];

FkCompute[configPath_String, opts : OptionsPattern[]] := Module[
  {pyExe, req},
  pyExe = FkComputePythonExe[OptionValue["PythonExecutable"]];
  req = <|"mode" -> "config", "config" -> configPath|>;
  FkComputeCallBridge[pyExe, req]
];

FkComputeVersion[opts : OptionsPattern[]] := Module[
  {pyExe, req},
  pyExe = FkComputePythonExe[OptionValue["PythonExecutable"]];
  req = <|"mode" -> "version"|>;
  FkComputeCallBridge[pyExe, req]
];

End[]
EndPackage[]

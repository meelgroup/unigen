{
  description = "Almost-uniform sampler";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    # cadiback, evalmaxsat and treedecomp are not used by unigen directly, but
    # our deps (cryptominisat, arjun, approxmc) pull them in, and approxmc's
    # exported config needs treedecomp. Declare them here as single anchors so
    # the whole tree shares one of each instead of duplicating them in flake.lock.
    cadiback = {
      url = "github:meelgroup/cadiback/main";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    evalmaxsat = {
      url = "github:meelgroup/EvalMaxSAT/master";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    treedecomp = {
      url = "github:meelgroup/treedecomp/main";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    cryptominisat = {
      url = "github:msoos/cryptominisat/master";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.cadiback.follows = "cadiback";
    };
    sbva = {
      url = "github:meelgroup/sbva/master";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    arjun = {
      url = "github:meelgroup/arjun/master";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.cadiback.follows = "cadiback";
      inputs.cadical.follows = "cadiback/cadical";
      inputs.cryptominisat.follows = "cryptominisat";
      inputs.sbva.follows = "sbva";
      inputs.evalmaxsat.follows = "evalmaxsat";
      inputs.treedecomp.follows = "treedecomp";
    };
    approxmc = {
      url = "github:meelgroup/approxmc/master";
      inputs.nixpkgs.follows = "nixpkgs";
      inputs.cadiback.follows = "cadiback";
      inputs.arjun.follows = "arjun";
      inputs.cryptominisat.follows = "cryptominisat";
      inputs.sbva.follows = "sbva";
      inputs.evalmaxsat.follows = "evalmaxsat";
      inputs.treedecomp.follows = "treedecomp";
    };
  };
  outputs =
    {
      self,
      nixpkgs,
      cadiback,
      evalmaxsat,
      treedecomp,
      cryptominisat,
      sbva,
      arjun,
      approxmc,
    }:
    let
      inherit (nixpkgs) lib;
      systems = lib.intersectLists lib.systems.flakeExposed lib.platforms.linux;
      forAllSystems = lib.genAttrs systems;
      nixpkgsFor = forAllSystems (system: nixpkgs.legacyPackages.${system});
      fs = lib.fileset;
      unigen-package =
        {
          stdenv,
          cmake,
          autoPatchelfHook,
          pkg-config,
          gmp,
          mpfr,
          zlib,
          cryptominisat,
          sbva,
          arjun,
          approxmc,
          evalmaxsat,
          treedecomp,
        }:
        stdenv.mkDerivation {
          name = "unigen";
          src = fs.toSource {
            root = ./.;
            fileset = fs.unions [
              ./CMakeLists.txt
              ./cmake
              ./src
              ./unigenConfig.cmake.in
            ];
          };
          nativeBuildInputs = [
            cmake
            autoPatchelfHook
            pkg-config
          ];
          cmakeFlags = [
            "-Dcryptominisat5_DIR=${cryptominisat}/lib/cmake/cryptominisat5"
            "-Darjun_DIR=${arjun}/lib/cmake/arjun"
            "-Dapproxmc_DIR=${approxmc}/lib/cmake/approxmc"
          ];
          buildInputs = [
            gmp
            mpfr
            zlib
            cryptominisat
            sbva
            arjun
            approxmc
            evalmaxsat
            treedecomp
          ];
        };
    in
    {
      packages = forAllSystems (
        system:
        let
          unigen = nixpkgsFor.${system}.callPackage unigen-package {
            cryptominisat = cryptominisat.packages.${system}.cryptominisat;
            sbva = sbva.packages.${system}.sbva;
            arjun = arjun.packages.${system}.arjun;
            approxmc = approxmc.packages.${system}.approxmc;
            evalmaxsat = evalmaxsat.packages.${system}.evalmaxsat;
            treedecomp = treedecomp.packages.${system}.treedecomp;
          };
        in
        {
          inherit unigen;
          default = unigen;
        }
      );
    };
}

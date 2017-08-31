/**********************************************************************
Copyright (C) 2017 by Noel M. O'Boyle

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/
#include "openbabel/babelconfig.h"
#include <string>
#include <algorithm>
#include "openbabel/mol.h"
#include "openbabel/obconversion.h"
#include "openbabel/reaction.h"

using namespace std;

namespace OpenBabel
{
  class ReactionInChIFormat : public OBFormat
  {
  public:
    //Register this format type ID
    ReactionInChIFormat()
    {
      OBConversion::RegisterFormat("rinchi",this);
    }

    virtual const char* Description()
    {
      return
        "Reaction SMILES format\n"
        "Write Options e.g. -xt\n"
        "  r radicals lower case eg ethyl is Cc\n"
        "\n";

    }

    virtual const char* GetMIMEType()
    { return "chemical/x-daylight-smiles"; }; // not right, need something else

    virtual const char* TargetClassDescription()
    {
      return OBReaction::ClassDescription();
    }

    const type_info& GetType()
    {
      return typeid(OBReaction*);
    }


    ////////////////////////////////////////////////////
    /// The "API" interface functions
    virtual bool ReadMolecule(OBBase* pReact, OBConversion* pConv);
    virtual bool WriteMolecule(OBBase* pReact, OBConversion* pConv);

    ////////////////////////////////////////////////////
    /// The "Convert" interface functions
    virtual bool ReadChemObject(OBConversion* pConv)
    {
      return true;
    }

    virtual bool WriteChemObject(OBConversion* pConv)
    {
      //WriteChemObject() always deletes the object retrieved by GetChemObject
      OBBase* pOb = pConv->GetChemObject();
      OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
      if(pReact==NULL)
        return false;

      bool ret=false;
      ret=WriteMolecule(pReact,pConv);

      delete pOb;
      return ret;
    }

  };

  //Make an instance of the format class
  ReactionInChIFormat theReactionInChIFormat;

  /////////////////////////////////////////////////////////////////
  bool ReactionInChIFormat::ReadMolecule(OBBase* pOb, OBConversion* pConv)
  {
    return true;
  }

  // Trim off the "InChI=1S/" and the trailing \n
  static std::string TrimInChI(const char *inchi)
  {
    std::string trimmed;
    const char *p = inchi + 9;
    while (true) {
      trimmed += *p;
      p++;
      if (*p == '\n' || *p == '\0')
        break;
    }
    return trimmed;
  }

  /////////////////////////////////////////////////////////////////
  bool ReactionInChIFormat::WriteMolecule(OBBase* pOb, OBConversion* pConv)
  {
    //It's really a reaction, not a molecule.
    //Cast output object to the class type need, i.e. OBReaction
    OBReaction* pReact = dynamic_cast<OBReaction*>(pOb);
    if(pReact==NULL)
      return false;
    ostream &ofs = *pConv->GetOutStream();

    OBFormat* pInChIFormat = OBConversion::FindFormat("inchi");
    if(!pInChIFormat)
      return false;

    OBConversion inchiconv;
    inchiconv.SetOutFormat(pInChIFormat);
    stringstream ss;
    inchiconv.SetOutStream(&ss);

#define REACTANTS 0
#define PRODUCTS 1
#define AGENTS 2

    std::vector<std::vector<std::string> > inchis(3);
    for (int part = 0; part <= 2; ++part) {
      unsigned int N;
      switch (part) {
      case REACTANTS: N = pReact->NumReactants(); break;
      case PRODUCTS: N = pReact->NumProducts(); break;
      case AGENTS: N = pReact->NumAgents(); break;
      }
      for (unsigned int i = 0; i < N; ++i) {
        OBMol* mol;
        switch (part) {
        case REACTANTS: mol = &*(pReact->GetReactant(i)); break;
        case PRODUCTS: mol = &*(pReact->GetProduct(i)); break;
        case AGENTS: mol = &*(pReact->GetAgent(i)); break;
        }
        bool ok = inchiconv.Write(mol);
        if (!ok)
          return false;

        string inchi = ss.str();
        if (strncmp(inchi.c_str(), "InChI=1S/", 9) != 0)
          return false;
        inchis[part].push_back(TrimInChI(inchi.c_str()));
        ss.str("");
      }
    }

    std::sort(inchis[REACTANTS].begin(), inchis[REACTANTS].end());
    std::sort(inchis[PRODUCTS].begin(), inchis[PRODUCTS].end());
    std::sort(inchis[AGENTS].begin(), inchis[AGENTS].end());

    bool reactants_first = true;
    std::vector<std::string>::const_iterator reactant_it = inchis[REACTANTS].begin();
    std::vector<std::string>::const_iterator product_it = inchis[PRODUCTS].end();
    while (true) {
      if (reactant_it == inchis[REACTANTS].end()) {
        if (product_it != inchis[PRODUCTS].end())
          reactants_first = false;
        break;
      }
      if (product_it == inchis[PRODUCTS].end()) 
        break;

      if (*reactant_it < *product_it) {
        reactants_first = false;
        break;
      }
      product_it++;
      reactant_it++;
    }

    ofs << "RInChI=1.00.1S/";
    std::vector<std::string>  &rxn_components = reactants_first ? inchis[REACTANTS] : inchis[PRODUCTS];
    for (std::vector<std::string>::const_iterator vit = rxn_components.begin(); vit != rxn_components.end(); ++vit) {
      if (vit != rxn_components.begin())
        ofs << '!';
      ofs << *vit;
    }
    ofs << "<>";
    rxn_components = reactants_first ? inchis[PRODUCTS] : inchis[REACTANTS];
    for (std::vector<std::string>::const_iterator vit = rxn_components.begin(); vit != rxn_components.end(); ++vit) {
      if (vit != rxn_components.begin())
        ofs << '!';
      ofs << *vit;
    }
    ofs << "<>";
    for (std::vector<std::string>::const_iterator vit = inchis[AGENTS].begin(); vit != inchis[AGENTS].end(); ++vit) {
      if (vit != inchis[AGENTS].begin())
        ofs << '!';
      ofs << *vit;
    }

    ofs << '\n';
    return true;
  }

} //namespace

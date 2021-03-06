

#include "MainFrame.hpp"

#include "Projectionist.hpp"
#include "FileManager.hpp"


#include <TGMsgBox.h>
#include <TGStatusBar.h>
#include <TGFileDialog.h>
#include <TGComboBox.h>
#include <TGTextEntry.h>

#include <TG3DLine.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGSlider.h>
#include <TH3D.h>
#include <TH3F.h>
#include <TH3I.h>

#include <TPython.h>
#include <Python.h>


#include <TString.h>

#include <array>
#include <tuple>
#include <vector>
#include <string>
#include <list>
#include <fstream>


ClassImp(MyMainFrame)



struct DataQuery {
  std::string cfg;
  std::string cent;
};

struct DataTree {

  std::vector<std::pair< std::string, // cfg
  std::vector<std::pair< std::string, // cent
  std::vector<std::pair< std::string, // kt
  std::vector<std::pair< std::string, // pair
  std::vector<std::pair< std::string, // magfield
                         std::string  // full-path
                         >> >> >> >> >> root;

  void insert_pieces(const std::vector<std::string> &l)
  {
    auto analysis = l[0],
         cfg = l[1],
         cent = l[2],
         kt = l[3],
         pair = l[4],
         magfield = l[5],
         path = l.back();

    path = analysis + "/" + cfg + "/" + pair + "/" + cent + "/" + kt + "/" + magfield;
    // tree.root[cfg][cent][kt][pair][magfield] = path;

    // auto found = std::find(root.begin(), root.end(), [cfg] (auto key) { return key->first == cfg; });
    std::pair<std::string, std::string> mp {magfield, path};

    std::pair<std::string, std::vector<decltype(mp)>> pairp {pair, {mp}};
    std::pair<std::string, std::vector<decltype(pairp)>> ktp {kt, {pairp}};
    std::pair<std::string, std::vector<decltype(ktp)>> centp {cent, {ktp}};
    std::pair<std::string, std::vector<decltype(centp)>> cfg_p {cfg, {centp}};

    // for (auto cfgit = root.begin(); cfgit != root.end(); ++cfgit) {
    for (auto &cfgit : root) {
      if (cfgit.first == cfg) {
      for (auto &cent_it : cfgit.second) {
        if (cent_it.first == cent) {
        for (auto &kt_it : cent_it.second) {
          if (kt_it.first == kt) {
          for (auto &pair_it : kt_it.second) {
            if (pair_it.first == pair) {
            for (auto &mfield_it : pair_it.second) {
              if (mfield_it.second == magfield) {
                std::cerr << "Unexpected duplicate path\n";
                return;
              }
            }
            pair_it.second.emplace_back(mp);
            return;
            }
          }
          kt_it.second.emplace_back(pairp);
          return;
          }
        }
        cent_it.second.emplace_back(ktp);
        return;
        }
      }
      cfgit.second.emplace_back(centp);
      return;
      }
    }
    root.emplace_back(cfg_p);
    return;
  }
};

struct FrameData {
  TString filename;
  TString sfname;

  TGLabel *lbl_fname;
  TGLabel *lbl_sname;

  TGLabel *lbl_xval,
          *lbl_yval,
          *lbl_zval;

  TGHSlider *slider_x,
            *slider_y,
            *slider_z;

  DataTree available_data;

  TGComboBox *cfg_ddown,
             *cent_ddown,
             *kt_ddown,
             *magfield_ddown,
             *pair_ddown;

  FileManager files;
  ProjectionistManager projection_manager;

  std::shared_ptr<TFile>
  open_file(const std::string filename)
    {
      auto tfile = std::make_shared<TFile>(filename.c_str());
      files.add(filename, tfile);
      return tfile;
    }

  void SetFitResultFilename(const TString &fname)
    {
      lbl_fname->SetText(fname.Data());
    }

  void SetSourceFilename(const TString &fname)
    {
      lbl_sname->SetText(fname.Data());
    }

  std::string get_selected_path()
    {
      const auto
        cfg_id = cfg_ddown->GetSelected(),
        cent_id = cent_ddown->GetSelected(),
        kt_id = kt_ddown->GetSelected(),
        pair_id = pair_ddown->GetSelected(),
        mfld_id = magfield_ddown->GetSelected();

      const auto cfgpair = available_data.root[cfg_id];
      const auto centpair = cfgpair.second[cent_id];
      const auto ktpair = centpair.second[kt_id];
      const auto pairpair = ktpair.second[pair_id];
      const auto mfield_pair = pairpair.second[mfld_id];

      return mfield_pair.second;
    }

  void build_slider_frame(MyMainFrame *self, TGCompositeFrame *parent)
    {
      TGHorizontalFrame *slidergroup = new TGHorizontalFrame(parent, 300, 50, kHorizontalFrame);
        slidergroup->SetName("slidergroup");

      auto add_labeled_slider = [&self](TGCompositeFrame *p, TString name)
        {
          TGVerticalFrame *frame = new TGVerticalFrame(p,100,50,kVerticalFrame);
            TGHorizontalFrame *lblframe = new TGHorizontalFrame(frame,318,19,kHorizontalFrame);
              TGLabel *axislbl = new TGLabel(lblframe, name);
                axislbl->SetName(Form("lbl%s", name.Data()));
                axislbl->SetTextJustify(36);
                axislbl->SetMargins(0,0,0,0);
                axislbl->SetWidth(100);
                axislbl->SetWrapLength(-1);
                lblframe->AddFrame(axislbl, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,0,2,0));

              TGLabel *lblval = new TGLabel(lblframe, "51");
                lblval->SetTextJustify(36);
                lblval->SetMargins(40,50,0,0);
                lblval->SetWrapLength(-1);
                lblframe->AddFrame(lblval, new TGLayoutHints(kLHintsLeft | kLHintsTop,2,2,2,0));
              frame->AddFrame(lblframe, new TGLayoutHints(kLHintsCenterX | kLHintsTop | kLHintsExpandX,2,2,2,0));

            TGHSlider *slider = new TGHSlider(frame, 0, kSlider1 | kScaleBoth, -1, kHorizontalFrame);
              slider->SetName(Form("slider_%s", name.Data()));
              slider->SetRange(0, 100);
              slider->SetPosition(50);
              frame->AddFrame(slider, new TGLayoutHints(kLHintsCenterX | kLHintsBottom | kLHintsExpandY| kLHintsExpandX,1,1,0,2));

            slider->Connect("PositionChanged(int)", "MyMainFrame", self, "OnSliderUpdate()");

           p->AddFrame(frame, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandX));
           return std::make_tuple(slider, lblval);
        };

      auto add_separator_line = [](TGCompositeFrame *p)
        {
          auto vline = new TGVertical3DLine(p);
          p->AddFrame(vline, new TGLayoutHints(kLHintsExpandY));
        };

      std::tie(slider_x, lbl_xval) = add_labeled_slider(slidergroup, "qOut");
      add_separator_line(slidergroup);
      std::tie(slider_y, lbl_yval) = add_labeled_slider(slidergroup, "qSide");
      add_separator_line(slidergroup);
      std::tie(slider_z, lbl_zval) = add_labeled_slider(slidergroup, "qLong");

      parent->AddFrame(slidergroup, new TGLayoutHints(kLHintsLeft | kLHintsExpandX | kLHintsTop,2,2,2,2));
      slidergroup->MoveResize(0,0,400,50);
    }

  void build_data_select_box(MyMainFrame *self, TGCompositeFrame *parent)
    {
      auto add_combo_box = [&](int id, TString label)
        {
          TGLabel *lbl = new TGLabel(parent, label);
            lbl->SetTextJustify(36);
            lbl->SetMargins(0,0,8,0);
            lbl->SetWrapLength(-1);
            parent->AddFrame(lbl, new TGLayoutHints(kLHintsLeft | kLHintsTop, 2,2,0,0));

          TGComboBox *cbox = new TGComboBox(parent, id, kHorizontalFrame | kSunkenFrame | kOwnBackground);
            cbox->Resize(247,25);
            parent->AddFrame(cbox, new TGLayoutHints(kLHintsNormal, 10));
            cbox->Connect("Selected(int, int)", "MyMainFrame", self, "OnDropdownSelection(int, int)");

          return cbox;
        };

      cfg_ddown = add_combo_box(1, "Config Hash");
      cent_ddown = add_combo_box(2, "Centrality");
      kt_ddown = add_combo_box(3, "kT");
      pair_ddown = add_combo_box(4, "Pair");
      magfield_ddown = add_combo_box(5, "Magnetic Field");
    }

  void UpdateLabels()
    {
      const Int_t
        Nx = slider_x->GetPosition(),
        Ny = slider_y->GetPosition(),
        Nz = slider_z->GetPosition();

      float frac_x, frac_y, frac_z;

      auto &r = projection_manager.current->ratio;
      if (r) {
        frac_x = r->GetXaxis()->GetBinCenter(Nx);
        frac_y = r->GetYaxis()->GetBinCenter(Ny);
        frac_z = r->GetZaxis()->GetBinCenter(Nz);
      } else {
        frac_x = float(slider_x->GetPosition()) / slider_x->GetMaxPosition();
        frac_y = float(slider_y->GetPosition()) / slider_y->GetMaxPosition();
        frac_z = float(slider_z->GetPosition()) / slider_z->GetMaxPosition();
      }

      lbl_xval->SetText(Form("%d  |  %.03f ", Nx, frac_x));
      lbl_yval->SetText(Form("%d  |  %.03f ", Ny, frac_y));
      lbl_zval->SetText(Form("%d  |  %.03f ", Nz, frac_z));
    }

  void update_choices(std::vector<std::string> cfg, std::vector<std::string> cent, std::vector<std::string> kt)
    {
      cent_ddown->RemoveAll();
      for (size_t i=0; i < cent.size(); ++i) {
        cent_ddown->AddEntry(cent[i].c_str(), i);
      }

      kt_ddown->RemoveAll();
      for (size_t i=0; i < kt.size(); ++i) {
        kt_ddown->AddEntry(kt[i].c_str(), i);
      }
    }

  void update_sliders()
    {
      auto update_slider = [](TGSlider &slider, TAxis &ax)
        {
          const Int_t
            N = ax.GetNbins(),
            M = slider.GetMaxPosition();

          if (N != M) {
            const Int_t
              n = slider.GetPosition(),
              val = N * double(n) / M + 1;

            slider.SetRange(1, N);
            slider.SetPosition(val);
          }
        };

      auto &r = projection_manager.current->ratio;
      if (!r) {
        return;
      }

      update_slider(*slider_x, *r->GetXaxis());
      update_slider(*slider_y, *r->GetYaxis());
      update_slider(*slider_z, *r->GetZaxis());

      UpdateLabels();
    }

  void reset_datatree(DataTree tree)
  {
    available_data = tree;

    cfg_ddown->RemoveAll();
    cent_ddown->RemoveAll();
    kt_ddown->RemoveAll();
    pair_ddown->RemoveAll();
    magfield_ddown->RemoveAll();

    int idx_cfg = 0;
    for (auto &cfg_pair : available_data.root) {
      cfg_ddown->AddEntry(cfg_pair.first.c_str(), idx_cfg++);

      if (idx_cfg == 1) {
        cfg_ddown->Select(0, false);

        int idx_cent = 0;
        for (auto &cent_pair : cfg_pair.second) {
          cent_ddown->AddEntry(cent_pair.first.c_str(), idx_cent++);
          if (idx_cent == 1) {
            cent_ddown->Select(0, false);

        int idx_kt = 0;
        for (auto &kt_pair : cent_pair.second) {
          kt_ddown->AddEntry(kt_pair.first.c_str(), idx_kt++);
          if (idx_kt == 1) {
            kt_ddown->Select(0, false);

        int idx_pair = 0;
        for (auto &pair_pair : kt_pair.second) {
          pair_ddown->AddEntry(pair_pair.first.c_str(), idx_pair++);
          if (idx_pair == 1) {
            pair_ddown->Select(0, false);

        int idx_magfield = 0;
        for (auto &magfield_pair : pair_pair.second) {
          magfield_ddown->AddEntry(magfield_pair.first.c_str(), idx_magfield++);
          if (idx_magfield == 1) {
            magfield_ddown->Select(0, true);
          } } } } } } }
        }
      }
    }
  }
};


MyMainFrame::MyMainFrame(const TGWindow *p)
  : TGMainFrame(p, 600, 400)
  , fCanvas(nullptr)
  , data(new FrameData())
{
  TGVerticalFrame *subframe = new TGVerticalFrame(this, 600, 400);

  TGHorizontalFrame *topbar = new TGHorizontalFrame(subframe, 200, 80);
    TGTextButton *open_button = new TGTextButton(topbar, "&Open");
    open_button->Connect("Clicked()", "MyMainFrame", this, "OnOpen()");
    topbar->AddFrame(open_button, new TGLayoutHints(kLHintsTop | kLHintsLeft, 0, 50));

    TGVerticalFrame *fname_frame = new TGVerticalFrame(topbar, 400, 80);
      data->lbl_fname = new TGLabel(fname_frame, "<nofile>");
      fname_frame->AddFrame(data->lbl_fname, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));
      data->lbl_fname->SetTextJustify(kTextLeft);

      data->lbl_sname = new TGLabel(fname_frame, "<nofile>");
      data->lbl_sname->SetTextJustify(kTextLeft);
      fname_frame->AddFrame(data->lbl_sname, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));
      topbar->AddFrame(fname_frame, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX));

    subframe->AddFrame(topbar, new TGLayoutHints(kLHintsExpandX));

  // { auto line = new TGHorizontal3DLine(topbar); topbar->AddFrame(line, new TGLayoutHints(kLHintsExpandX)); }

  TGHorizontalFrame *mainframe = new TGHorizontalFrame(subframe, 600, 200);
    TGVerticalFrame *leftbar = new TGVerticalFrame(mainframe);
      TGGroupFrame *databox = new TGGroupFrame(leftbar, "Data Select");
      data->build_data_select_box(this, databox);
      leftbar->AddFrame(databox);

      fitresult_panel = std::make_unique<PanelFitResult>(*leftbar, kLHintsLeft | kLHintsTop | kLHintsExpandX);
      cutconfig_panel = std::make_unique<PanelCutConfig>(*leftbar);

      mainframe->AddFrame(leftbar, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY, 2,2,2,2));

    auto *slider_canvas_container = new TGVerticalFrame(mainframe, 400, 200);
      data->build_slider_frame(this, slider_canvas_container);

      fCanvas = new TRootEmbeddedCanvas(0, slider_canvas_container, 500, 400);
      Int_t wid = fCanvas->GetCanvasWindowId();
      auto canvas = new TCanvas("MyCanvas", 10, 10, wid);
      // canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","MyMainFrame",this,
      //              "EventInfo(Int_t,Int_t,Int_t,TObject*)");

      fCanvas->AdoptCanvas(canvas);

      auto canvas_props = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX | kLHintsExpandY,0,0,1,1);
      slider_canvas_container->AddFrame(fCanvas, canvas_props);

    mainframe->AddFrame(slider_canvas_container, new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsExpandY | kLHintsExpandX, 2,2,2,2));

  subframe->AddFrame(mainframe, new TGLayoutHints(kLHintsExpandX  | kLHintsExpandY, 0, 0, 0, 0));
  Int_t parts[] = {45, 15, 10, 30};
  auto fStatusBar = new TGStatusBar(subframe, 50, 10, kVerticalFrame);
  fStatusBar->SetParts(parts, 4);
  fStatusBar->Draw3DCorner(kFALSE);
  subframe->AddFrame(fStatusBar, new TGLayoutHints(kLHintsExpandX, 0, 0, 10, 0));

  TGHorizontalFrame *hframe = new TGHorizontalFrame(subframe, 10, 40);

   TGTextButton *draw = new TGTextButton(hframe, "&Draw");
   draw->Connect("Clicked()", "MyMainFrame", this, "OnDraw()");

   hframe->AddFrame(draw, new TGLayoutHints(kLHintsCenterX, 5, 5, 3, 4));
   subframe->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX, 2, 2, 2, 2));


  //  AddFrame(draw, new TGTextButton(hframe, "&Draw "))

  AddFrame(subframe, new TGLayoutHints(kLHintsCenterX | kLHintsCenterY | kLHintsExpandX | kLHintsExpandY));

  SetWindowName("FitInspector");

  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
  Resize(1266, 1051);
}


MyMainFrame::~MyMainFrame()
{
  Cleanup();
  delete fCanvas;
}

void
MyMainFrame::OnDraw()
{
  std::cout << "OnDraw\n";

  TCanvas &c = *fCanvas->GetCanvas();

  c.SetFillColor(42);

  c.Modified();
  c.Update();
}


void
MyMainFrame::OnOpen()
{
  std::cout << "OnOpen\n";

  TGFileInfo fi;
  static const char* filetypes[] = {
    "JSON files", "*.json",
    "All files", "*",
    nullptr, nullptr
    };

  fi.fFileTypes = filetypes;

  new TGFileDialog(
    gClient->GetRoot(),
    this,
    kFDOpen,
    &fi);

  TString selected_filename(fi.fFilename);

  // no selection
  if (selected_filename.IsWhitespace()) {
    return;
  }

  auto result = LoadJsonFile(selected_filename);
  if (result.is_error) {
    std::cerr << result.err.reason << "\n";
    new TGMsgBox(gClient->GetRoot(),
                 gClient->GetRoot(),
                 "Error",
                 result.err.reason.Data(),
                 EMsgBoxIcon::kMBIconExclamation);
  }
}

JsonLoadResult
MyMainFrame::LoadJsonFile(TString filename)
{
  auto pyobj_to_str = [] (PyObject *obj)
    {
      Py_ssize_t len = 0;
      const char *s = PyUnicode_AsUTF8AndSize(obj, &len);
      return std::string(s, s+len);
    };

  auto pyobj_to_vec_str = [pyobj_to_str] (PyObject *pylist)
    {
      std::vector<std::string> result;
      auto count = PyList_Size(pylist);
      for (int i=0; i<count; ++i) {
        auto pystr = PyList_GetItem(pylist, i);
        result.emplace_back(pyobj_to_str(pystr));
      }
      return result;
    };

  auto load_string_vec = [pyobj_to_vec_str] (TString cmd)
    {
      PyObject *pylist = TPython::Eval(cmd);
      std::vector<std::string> result = pyobj_to_vec_str(pylist);
      Py_DECREF(pylist);
      return result;
    };

  auto load_string_vec_vec = [pyobj_to_vec_str] (TString cmd)
    {
      std::vector<std::vector<std::string>> result;
      PyObject *pylist = TPython::Eval(cmd);
      auto count = PyList_Size(pylist);
      for (int i=0; i<count; ++i) {
        auto pylist_inner = PyList_GetItem(pylist, i);

        std::vector<std::string> inresult = pyobj_to_vec_str(pylist_inner);
        result.emplace_back(std::move(inresult));

        Py_DECREF(pylist_inner);
      }
      Py_DECREF(pylist);
      return result;
    };

  // try open-file
  std::ifstream file(filename);
  if (!file.is_open()) {
    return JsonLoadResult::Error(Form("Could not open file '%s'", filename.Data()));
  }

  if (!TPython::Exec("import json; from pathlib import Path")) {
    return JsonLoadResult::Error("Could not import modules from standard libray");
  }

  if (!TPython::Exec("import pandas as pd")) {
    return JsonLoadResult::Error("Could not import pandas");
  }

  const char *cmd = Form("data = json.loads(Path('%s').read_text())", filename.Data());
  if (!TPython::Exec(cmd)) {
    return JsonLoadResult::Error("Error Reading " + filename + " -- Expected valid JSON file");
  }

  const PyObject *data_res = TPython::Eval("data");
  if (data_res == nullptr) {
    return JsonLoadResult::Error("Error Reading " + filename + " -- Expected valid JSON file");
  }

  auto *df_res = PyDict_GetItemString(const_cast<PyObject*>(data_res), "df");
  if (df_res == nullptr) {
    return JsonLoadResult::Error("Invalid file " + filename + " -- missing 'df' key");
  }

  auto *filename_res = PyDict_GetItemString(const_cast<PyObject*>(data_res), "filename");
  if (filename_res == nullptr) {
    return JsonLoadResult::Error("Invalid file " + filename + " -- missing 'filename' key");
  }

  if (!PyUnicode_Check(filename_res)) {
    return JsonLoadResult::Error("Path to rootfile not a string");
  }

  std::string root_file = pyobj_to_str(filename_res);

  std::shared_ptr<TFile> rfile = data->open_file(root_file);
  if (rfile == nullptr) {
    return JsonLoadResult::Error("Could not read data-file " + root_file + " -- is it a valid ROOT file?");
  }

  if (!TPython::Exec("df = pd.DataFrame(data['df'])")) {
    return JsonLoadResult::Error("Could not read dataframe in file " + root_file);
  }

  // remove "bad" entries
  if (!TPython::Exec("df = df[~df.kt.isnull()]")) {
    return JsonLoadResult::Error("Bad dataframe in file " + root_file);
  }

  data->SetSourceFilename(root_file);
  data->SetFitResultFilename(filename.Data());

  std::vector<std::string> cent = load_string_vec("list(df.cent.unique())");
  std::vector<std::string> kt = load_string_vec("list(df.kt.unique())");

  DataTree tree;

  auto vv = load_string_vec_vec("[(list(l) + ['/'.join(l)])"
                                " for l in df[['analysis', 'cfg', 'cent', 'kt', 'pair', 'magfield']].drop_duplicates().values]");

  for (auto &l : vv) {
    tree.insert_pieces(l);
  }

  data->reset_datatree(tree);

  return JsonLoadResult::Ok(Form("Loaded file '%s'", filename.Data()));
}

void
MyMainFrame::OnSliderUpdate()
{
  const std::array<Int_t,3>
    idxs = {data->slider_x->GetPosition(),
            data->slider_y->GetPosition(),
            data->slider_z->GetPosition()};

  auto projectionist = data->projection_manager.current;
  auto pad = projectionist->get_full_canvas(idxs);

  TCanvas &canvas = *fCanvas->GetCanvas();
  canvas.Clear();
  canvas.cd();

  pad->ResizePad();
  pad->Draw();
  canvas.Modified();
  canvas.Update();

  data->UpdateLabels();
}

void
MyMainFrame::OnDropdownSelection(int id, int entry)
{
  auto path = data->get_selected_path();
  std::cout << " loading " << path << "\n";

  auto file = data->files.current_file;

  auto *tdir = dynamic_cast<TDirectory*>(file->Get(path.c_str()));
  if (!tdir) {
    std::cerr << "ERROR: Could not load path '" << path << "'\n";
    return;
  }

  if (id == 1) {
    std::cout << "Updating cut configuration\n";
  } else {
    std::cout << "Updating data selection\n";
  }

  auto &tree = data->available_data.root;
  auto &cfg_node = tree[data->cfg_ddown->GetSelected()];
  auto &cent_node = cfg_node.second[data->cent_ddown->GetSelected()];
  auto &kt_node = cent_node.second[data->kt_ddown->GetSelected()];
  auto &pair_node = kt_node.second[data->pair_ddown->GetSelected()];
  auto &magfield_node = pair_node.second[data->magfield_ddown->GetSelected()];

  const std::string
    &cfg = cfg_node.first,
    &cent = cent_node.first,
    &kt = kt_node.first,
    &pair = pair_node.first,
    &mag = magfield_node.first;

  auto *cmd = Form("series = df[(df.cfg=='%s')&(df.cent=='%s')&(df.kt=='%s')&(df.pair=='%s')&(df.magfield=='%s')&df.mrc.isnull()].iloc[0]",
                   cfg.c_str(),
                   cent.c_str(),
                   kt.c_str(),
                   pair.c_str(),
                   mag.c_str());

  // std::cout << "> " << cmd << "\n";
  if (!TPython::Exec(cmd)) {
    std::cerr << "ERROR: Could not load fit results from selected data" << cmd << "\n";
  } else {
    auto *series = static_cast<PyObject*>(TPython::Eval("dict(series)"));
    std::cout << "Loaded series " << series << "\n";

    std::vector<std::pair<std::string, double>> values;
    std::vector<const char*> keys = {"Ro", "Rs", "Rl", "lam", "norm"};
    for (auto &key : keys) {
      auto *pyval = PyDict_GetItemString(series, key);
      if (PyFloat_Check(pyval)) {
        double value = PyFloat_AS_DOUBLE(pyval);
        values.emplace_back(key, value);
      }
    }

    fitresult_panel->Update(values);
  }

  std::cout << "ID: " << id << ", " << entry << "\n";

  // update configuration panel
  if (id == 5 || id == 1) {
    const TString cfg_path = "Q3DLCMS/" + cfg;
    if (auto *dir = dynamic_cast<TDirectory*>(file->Get(cfg_path))) {
      cutconfig_panel->Update(*dir);
    }
    else {
      std::cout << "Could not find directory " << cfg_path << "\n";
    }
  }

  data->projection_manager.add_tdir(path, *tdir);
  data->update_sliders();
  OnSliderUpdate();
}
